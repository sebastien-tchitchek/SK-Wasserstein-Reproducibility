#!/usr/bin/env python3


from __future__ import annotations

import argparse
import http.cookiejar
import os
import shutil
import sys
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, Optional, Tuple


class DownloadFormParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.current_action: Optional[str] = None
        self.current_inputs: Dict[str, str] = {}
        self.best: Optional[Tuple[str, Dict[str, str]]] = None

    def handle_starttag(self, tag: str, attrs) -> None:
        data = dict(attrs)
        if tag.lower() == "form":
            self.current_action = data.get("action")
            self.current_inputs = {}
        elif tag.lower() == "input" and self.current_action:
            name = data.get("name")
            if name:
                self.current_inputs[name] = data.get("value", "")

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() != "form" or not self.current_action:
            return
        if "id" in self.current_inputs or "uuid" in self.current_inputs:
            self.best = (self.current_action, dict(self.current_inputs))
        self.current_action = None
        self.current_inputs = {}


def open_request(opener, url: str):
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                "Chrome/124 Safari/537.36"
            )
        },
    )
    return opener.open(request, timeout=120)


def response_is_file(response) -> bool:
    content_disposition = response.headers.get("Content-Disposition", "")
    content_type = response.headers.get_content_type()
    return "attachment" in content_disposition.lower() or content_type not in {
        "text/html",
        "text/plain",
    }


def stream_response(response, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    with temporary.open("wb") as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)
    os.replace(temporary, destination)


def download(file_id: str, destination: Path) -> None:
    cookie_jar = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(
        urllib.request.HTTPCookieProcessor(cookie_jar)
    )
    urls = (
        "https://drive.usercontent.google.com/download?"
        + urllib.parse.urlencode(
            {"id": file_id, "export": "download", "confirm": "t"}
        ),
        "https://docs.google.com/uc?"
        + urllib.parse.urlencode(
            {"id": file_id, "export": "download", "confirm": "t"}
        ),
    )

    last_error: Optional[Exception] = None
    for initial_url in urls:
        try:
            response = open_request(opener, initial_url)
            if response_is_file(response):
                stream_response(response, destination)
                return

            html = response.read(5 * 1024 * 1024).decode("utf-8", errors="replace")
            parser = DownloadFormParser()
            parser.feed(html)
            if parser.best is None:
                raise RuntimeError(
                    "Google Drive did not provide the download form."
                )
            action, fields = parser.best
            action = urllib.parse.urljoin(initial_url, action)
            final_url = action + "?" + urllib.parse.urlencode(fields)
            response = open_request(opener, final_url)
            if not response_is_file(response):
                text = response.read(2000).decode("utf-8", errors="replace")
                raise RuntimeError(
                    "Unexpected Google Drive response: " + text[:300].replace("\n", " ")
                )
            stream_response(response, destination)
            return
        except Exception as exc:  
            last_error = exc

    raise RuntimeError(f"Download failed: {last_error}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("file_id")
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    try:
        download(args.file_id, args.destination.expanduser().resolve())
    except Exception as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
