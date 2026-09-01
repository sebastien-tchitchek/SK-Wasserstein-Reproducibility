#!/usr/bin/env python3


from __future__ import annotations

import argparse
import bz2
import errno
import gzip
import hashlib
import lzma
import os
import shutil
import subprocess
import sys
import tarfile
import traceback
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

VTK_SUFFIXES = {
    ".vtk", ".vti", ".vtu", ".vtp", ".vtr", ".vts", ".vtm",
    ".pvti", ".pvtu", ".pvtp", ".pvtr", ".pvts",
}


MAGIC_XZ = b"\xfd7zXZ\x00"
MAGIC_GZIP = b"\x1f\x8b"
MAGIC_BZIP2 = b"BZh"
MAGIC_ZIP_PREFIXES = (b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")
MAGIC_7Z = b"7z\xbc\xaf'\x1c"
MAGIC_ZSTD = b"\x28\xb5\x2f\xfd"


@dataclass(frozen=True)
class ArchiveInfo:
    path: Path
    kind: str
    size: int


def human_size(value: int) -> str:
    value_f = float(max(0, value))
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    index = 0
    while value_f >= 1024.0 and index < len(units) - 1:
        value_f /= 1024.0
        index += 1
    return f"{value_f:.2f} {units[index]}"


def read_magic(path: Path, count: int = 512) -> bytes:
    try:
        with path.open("rb") as stream:
            return stream.read(count)
    except OSError:
        return b""


def looks_like_uncompressed_tar(header: bytes) -> bool:
    return len(header) >= 265 and header[257:262] in {b"ustar", b"ustar\x00"}


def archive_kind(path: Path) -> Optional[str]:
    
    try:
        if not path.is_file() or path.is_symlink() or path.stat().st_size < 2:
            return None
    except OSError:
        return None

    suffix = path.suffix.casefold()
    if suffix in VTK_SUFFIXES:
        return None
    if path.name.startswith(".") and path.stat().st_size < 1024 * 1024:
        return None

    
    
    try:
        if zipfile.is_zipfile(path):
            return "zip"
    except OSError:
        pass
    try:
        if tarfile.is_tarfile(path):
            return "tar"
    except (OSError, tarfile.TarError, EOFError, lzma.LZMAError):
        pass

    header = read_magic(path)
    if header.startswith(MAGIC_XZ):
        return "xz-single"
    if header.startswith(MAGIC_GZIP):
        return "gzip-single"
    if header.startswith(MAGIC_BZIP2):
        return "bzip2-single"
    if any(header.startswith(prefix) for prefix in MAGIC_ZIP_PREFIXES):
        return "zip"
    if header.startswith(MAGIC_7Z):
        return "7z"
    if header.startswith(MAGIC_ZSTD):
        return "zstd"
    if looks_like_uncompressed_tar(header):
        return "tar"
    return None


def iter_regular_files(root: Path) -> list[Path]:
    result: list[Path] = []
    for current, dirs, files in os.walk(root, followlinks=False):
        dirs[:] = sorted(d for d in dirs if not d.startswith(".wmt_partiel_"))
        base = Path(current)
        for name in sorted(files):
            path = base / name
            try:
                if path.is_file() and not path.is_symlink():
                    result.append(path)
            except OSError:
                continue
    return result


def find_archives(root: Path) -> list[ArchiveInfo]:
    result: list[ArchiveInfo] = []
    for path in iter_regular_files(root):
        kind = archive_kind(path)
        if kind is None:
            continue
        try:
            size = path.stat().st_size
        except OSError:
            size = 0
        result.append(ArchiveInfo(path=path, kind=kind, size=size))
    result.sort(key=lambda item: (len(item.path.parts), str(item.path).casefold()))
    return result


def find_vtk_files(root: Path) -> list[Path]:
    return [
        path for path in iter_regular_files(root)
        if path.suffix.casefold() in VTK_SUFFIXES
    ]


def xz_uncompressed_size(path: Path) -> Optional[int]:
    if not read_magic(path, 8).startswith(MAGIC_XZ):
        return None
    try:
        completed = subprocess.run(
            ["xz", "--robot", "--list", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    for line in reversed(completed.stdout.splitlines()):
        fields = line.split("\t")
        if fields and fields[0] == "totals" and len(fields) >= 5:
            try:
                return int(fields[4])
            except ValueError:
                return None
    return None


def zip_uncompressed_size(path: Path) -> Optional[int]:
    try:
        with zipfile.ZipFile(path) as archive:
            return sum(info.file_size for info in archive.infolist())
    except (OSError, zipfile.BadZipFile):
        return None


def estimated_uncompressed_size(info: ArchiveInfo) -> Optional[int]:
    if info.kind == "zip":
        return zip_uncompressed_size(info.path)
    if info.kind in {"tar", "xz-single"}:
        estimate = xz_uncompressed_size(info.path)
        if estimate is not None:
            return estimate
    if info.kind == "tar" and read_magic(info.path, 512)[257:262] in {b"ustar", b"ustar\x00"}:
        return info.size
    return None


def ensure_disk_space(info: ArchiveInfo) -> None:
    estimate = estimated_uncompressed_size(info)
    if estimate is None:
        return
    free = shutil.disk_usage(info.path.parent).free
    reserve = 1024**3  
    if free < estimate + reserve:
        raise RuntimeError(
            "Insufficient disk space to extract "
            f"{info.path.name}. Libre : {human_size(free)} ; "
            f"estimated requirement including margin: {human_size(estimate + reserve)}."
        )


def safe_target_for(path: Path) -> Path:
    name = path.name
    for ending in (
        ".tar.xz", ".tar.gz", ".tar.bz2", ".tar.zst",
        ".txz", ".tgz", ".tbz2", ".zip", ".7z", ".tar",
        ".xz", ".gz", ".bz2", ".zst",
    ):
        if name.casefold().endswith(ending):
            name = name[: -len(ending)]
            break
    cleaned = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)
    cleaned = cleaned.strip("._-") or "archive"
    digest = hashlib.sha1(str(path.absolute()).encode("utf-8")).hexdigest()[:8]
    return path.parent / f"{cleaned}__extracted_{digest}"


def _safe_member_destination(target: Path, member_name: str) -> Path:
    candidate = (target / member_name).resolve()
    target_resolved = target.resolve()
    try:
        candidate.relative_to(target_resolved)
    except ValueError as exc:
        raise RuntimeError(
            f"Unsafe path in archive: {member_name!r}"
        ) from exc
    return candidate


def extract_tar(path: Path, target: Path) -> None:
    with tarfile.open(path, mode="r:*") as archive:
        
        
        try:
            archive.extractall(target, filter="data")
        except TypeError:
            members = archive.getmembers()
            for member in members:
                _safe_member_destination(target, member.name)
                if member.issym() or member.islnk():
                    
                    
                    
                    continue
                archive.extract(member, target)


def extract_zip(path: Path, target: Path) -> None:
    with zipfile.ZipFile(path) as archive:
        for member in archive.infolist():
            _safe_member_destination(target, member.filename)
        archive.extractall(target)


def copy_decompressed(source, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as output:
        shutil.copyfileobj(source, output, length=16 * 1024 * 1024)


def extract_single_compressed(path: Path, target: Path, kind: str) -> None:
    target.mkdir(parents=True, exist_ok=True)
    name = path.name
    suffix_map = {
        "xz-single": ".xz",
        "gzip-single": ".gz",
        "bzip2-single": ".bz2",
    }
    suffix = suffix_map[kind]
    if name.casefold().endswith(suffix):
        output_name = name[: -len(suffix)] or "decompressed_content"
    else:
        output_name = name + ".decompresse"
    destination = target / output_name
    opener = {
        "xz-single": lzma.open,
        "gzip-single": gzip.open,
        "bzip2-single": bz2.open,
    }[kind]
    with opener(path, "rb") as source:
        copy_decompressed(source, destination)


def extract_external(path: Path, target: Path, kind: str) -> None:
    target.mkdir(parents=True, exist_ok=True)
    if kind == "7z":
        command = shutil.which("7z") or shutil.which("7zz")
        if command is None:
            raise RuntimeError(
                "A 7z archive was found, but the 7z tool is not installed. "
                "Install p7zip-full, then rerun bash run.sh."
            )
        subprocess.run([command, "x", "-y", f"-o{target}", str(path)], check=True)
        return
    if kind == "zstd":
        
        if shutil.which("zstd") is None:
            raise RuntimeError(
                "A zstd archive was found, but zstd is not installed. "
                "Install zstd, then rerun bash run.sh."
            )
        result = subprocess.run(
            ["tar", "--zstd", "-xf", str(path), "-C", str(target)],
            check=False,
        )
        if result.returncode == 0:
            return
        destination = target / (path.stem or "decompressed_content")
        with destination.open("wb") as output:
            subprocess.run(["zstd", "-dc", str(path)], check=True, stdout=output)
        return
    raise RuntimeError(f"Unsupported archive type: {kind}")


def extract_one(info: ArchiveInfo) -> tuple[Path, int]:
    path = info.path
    target = safe_target_for(path)
    marker = target / ".WMT_EXTRACTION_COMPLETE"

    if marker.is_file():
        
        
        try:
            path.unlink()
        except OSError:
            pass
        return target, len(iter_regular_files(target))

    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True, exist_ok=True)
    ensure_disk_space(info)

    try:
        if info.kind == "tar":
            extract_tar(path, target)
        elif info.kind == "zip":
            extract_zip(path, target)
        elif info.kind in {"xz-single", "gzip-single", "bzip2-single"}:
            extract_single_compressed(path, target, info.kind)
        elif info.kind in {"7z", "zstd"}:
            extract_external(path, target, info.kind)
        else:
            raise RuntimeError(f"Unknown type: {info.kind}")

        files = [p for p in iter_regular_files(target) if p.name != marker.name]
        if not files:
            raise RuntimeError("the archive produced no file")
        marker.write_text(
            f"Source archive: {path}\nExtracted files: {len(files)}\n",
            encoding="utf-8",
        )
        
        path.unlink()
        return target, len(files)
    except OSError as exc:
        if exc.errno == errno.ENOSPC:
            raise RuntimeError(
                "The disk became full during extraction. Free some "
                "space, then rerun exactly bash run.sh."
            ) from exc
        raise
    except Exception:
        shutil.rmtree(target, ignore_errors=True)
        raise


def write_inventory(report, root: Path, title: str) -> None:
    files = iter_regular_files(root)
    report.write(f"\n{title}\n")
    report.write(f"Number of files: {len(files)}\n")
    for path in files[:500]:
        try:
            size = path.stat().st_size
        except OSError:
            size = 0
        kind = archive_kind(path) or (
            "vtk" if path.suffix.casefold() in VTK_SUFFIXES else "file"
        )
        try:
            relative = path.relative_to(root)
        except ValueError:
            relative = path
        report.write(f"- [{kind}] {human_size(size):>12}  {relative}\n")
    if len(files) > 500:
        report.write(f"... {len(files) - 500} additional files not shown\n")


def prepare(
    root: Path,
    outer_archive: Optional[Path],
    report_path: Path,
) -> int:
    root = root.expanduser().resolve()
    report_path = report_path.expanduser().resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    marker = root / ".nested_archives_prepared"

    with report_path.open("w", encoding="utf-8") as report:
        report.write("PREPARATION OF NESTED WASSERSTEINMERGETREESDATA ARCHIVES\n")
        report.write(f"Root: {root}\n")
        write_inventory(report, root, "INITIAL INVENTORY")
        report.flush()

        initial_archives = find_archives(root)
        initial_vtk = find_vtk_files(root)
        print(
            f"  Current layer: {len(initial_archives)} archive(s), "
            f"{len(initial_vtk)} VTK file(s)."
        )

        
        
        
        if outer_archive is not None:
            outer_archive = outer_archive.expanduser().resolve()
            if outer_archive.is_file() and (initial_archives or initial_vtk):
                size = outer_archive.stat().st_size
                outer_archive.unlink()
                print(
                    "  Outer archive no longer needed and deleted "
                    f"({human_size(size)} freed)."
                )
                report.write(
                    f"\nOuter archive deleted: {outer_archive} "
                    f"({human_size(size)})\n"
                )

        extracted_total = 0
        for pass_number in range(1, 33):
            archives = find_archives(root)
            if not archives:
                break
            print(
                f"  Extraction level {pass_number} : "
                f"{len(archives)} archive(s) to open."
            )
            report.write(
                f"\nLEVEL {pass_number} — {len(archives)} ARCHIVE(S)\n"
            )
            report.flush()

            for index, info in enumerate(archives, start=1):
                try:
                    relative = info.path.relative_to(root)
                except ValueError:
                    relative = info.path
                print(
                    f"    [{index}/{len(archives)}] {relative} "
                    f"({info.kind}, {human_size(info.size)})"
                )
                report.write(
                    f"- extraction [{info.kind}] {human_size(info.size)}: "
                    f"{relative}\n"
                )
                report.flush()
                target, count = extract_one(info)
                extracted_total += 1
                report.write(
                    f"  -> {target.relative_to(root)} ({count} files)\n"
                )
                report.flush()
        else:
            raise RuntimeError(
                "Too many nested archive levels (more than 32)."
            )

        remaining = find_archives(root)
        vtk_files = find_vtk_files(root)
        write_inventory(report, root, "FINAL INVENTORY")
        report.write(f"\nArchives extracted successfully: {extracted_total}\n")
        report.write(f"Remaining recognized archives: {len(remaining)}\n")
        report.write(f"VTK files found: {len(vtk_files)}\n")

        if remaining:
            names = ", ".join(str(item.path.relative_to(root)) for item in remaining[:20])
            raise RuntimeError(
                "Some recognized archives were not extracted: " + names
            )
        if not vtk_files:
            raise RuntimeError(
                "No VTK file was found after extracting "
                "all layers. See the report: " + str(report_path)
            )

        marker.write_text(
            f"VTK files found: {len(vtk_files)}\n"
            f"Extracted archives: {extracted_total}\n",
            encoding="utf-8",
        )
        print(
            f"  Preparation completed: {len(vtk_files)} VTK file(s) available."
        )
        print(f"  Complete inventory: {report_path}")
        return len(vtk_files)


def main(argv: Sequence[str] = sys.argv[1:]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--outer-archive", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args(argv)

    try:
        prepare(args.data_root, args.outer_archive, args.report)
        return 0
    except Exception as exc:
        print(f"INTERNAL EXTRACTION ERROR: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
