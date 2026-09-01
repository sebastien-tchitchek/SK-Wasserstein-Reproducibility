#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PARAVIEW_REPOSITORY="https://github.com/topology-tool-kit/ttk-paraview.git"
PARAVIEW_TAG="v5.13.0"
PARAVIEW_COMMIT="8b383eeb53821e71f08593ac431a1b3b1855ccac"
PARAVIEW_SOURCE="$ROOT/third_party/ttk-paraview"
TTK_SOURCE="$ROOT/third_party/ttk"
PARAVIEW_BUILD="$ROOT/build/ttk-paraview"
TTK_BUILD="$ROOT/build/ttk"
PARAVIEW_JOBS="${SKW_BUILD_JOBS:-2}"
TTK_JOBS="${SKW_TTK_BUILD_JOBS:-1}"

TTK_DEB_PACKAGES="$(dpkg -l | awk '/^ii/ && ($2 == "ttk" || $2 == "ttk-paraview" || $2 == "paraview" || $2 == "python3-paraview") {print $2}')"
if [ -n "$TTK_DEB_PACKAGES" ]; then
  sudo apt-get remove --purge -y $TTK_DEB_PACKAGES
fi

sudo apt-get update
sudo apt-get install -y software-properties-common
sudo add-apt-repository -y universe
sudo apt-get update
sudo apt-get install -y \
  build-essential ca-certificates git cmake ninja-build pkg-config \
  libboost-system-dev libopengl-dev libgl1-mesa-dev libxcursor-dev libxt-dev \
  libx11-dev libxext-dev libxrender-dev libxfixes-dev libxrandr-dev \
  libxinerama-dev libxi-dev libxkbcommon-x11-dev \
  qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools qttools5-dev \
  qtxmlpatterns5-dev-tools libqt5xmlpatterns5-dev qtdeclarative5-dev \
  libqt5x11extras5-dev libqt5svg5-dev \
  libeigen3-dev libgraphviz-dev graphviz libqhull-dev \
  python3 python3-dev python3-numpy python3-sklearn \
  libsqlite3-dev libwebsocketpp-dev zlib1g-dev

rm -rf "$PARAVIEW_SOURCE" "$PARAVIEW_BUILD" "$TTK_BUILD"
mkdir -p "$ROOT/build"

git clone \
  --depth 1 \
  --branch "$PARAVIEW_TAG" \
  --single-branch \
  --recurse-submodules \
  --shallow-submodules \
  "$PARAVIEW_REPOSITORY" \
  "$PARAVIEW_SOURCE"

if [ "$(git -C "$PARAVIEW_SOURCE" rev-parse HEAD)" != "$PARAVIEW_COMMIT" ]; then
  echo "Unexpected ttk-paraview commit." >&2
  return 1 2>/dev/null || exit 1
fi

cmake -S "$PARAVIEW_SOURCE" -B "$PARAVIEW_BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr/local \
  -DPARAVIEW_PYTHON_SITE_PACKAGES_SUFFIX=lib/python3.12/site-packages \
  -DPARAVIEW_USE_PYTHON=ON \
  -DPARAVIEW_USE_QT=ON \
  -DPARAVIEW_USE_MPI=OFF \
  -DBUILD_TESTING=OFF

cmake --build "$PARAVIEW_BUILD" --parallel "$PARAVIEW_JOBS"
sudo rm -f \
  /usr/local/bin/paraview \
  /usr/local/bin/pvpython \
  /usr/local/bin/pvbatch \
  /usr/local/bin/pvserver
sudo cmake --install "$PARAVIEW_BUILD"

cmake -S "$TTK_SOURCE" -B "$TTK_BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr/local \
  -DCMAKE_PREFIX_PATH=/usr/local \
  -DParaView_DIR=/usr/local/lib/cmake/paraview-5.13 \
  -DTTK_INSTALL_PLUGIN_DIR=bin/plugins \
  -DTTK_PYTHON_MODULE_DIR=lib/python3.12/site-packages \
  -DTTK_BUILD_PARAVIEW_PLUGINS=ON \
  -DTTK_BUILD_VTK_WRAPPERS=ON \
  -DTTK_BUILD_VTK_PYTHON_MODULE=ON \
  -DTTK_BUILD_STANDALONE_APPS=ON \
  -DTTK_BUILD_DOCUMENTATION=OFF \
  -DTTK_ENABLE_KAMIKAZE=ON \
  -DTTK_ENABLE_DOUBLE_TEMPLATING=OFF \
  -DTTK_ENABLE_CPU_OPTIMIZATION=OFF \
  -DTTK_ENABLE_SHARED_BASE_LIBRARIES=ON \
  -DTTK_ENABLE_OPENMP=ON \
  -DTTK_ENABLE_MPI=OFF

cmake --build "$TTK_BUILD" --parallel "$TTK_JOBS"
sudo rm -rf /usr/local/bin/plugins/TopologyToolKit
sudo cmake --install "$TTK_BUILD"
sudo ldconfig

BASHRC="$HOME/.bashrc"
touch "$BASHRC"
for LINE in \
  'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib' \
  'export PV_PLUGIN_PATH=${PV_PLUGIN_PATH}:/usr/local/bin/plugins/' \
  'export PYTHONPATH=${PYTHONPATH}:/usr/local/lib/python3.12/site-packages/'
do
  grep -qxF "$LINE" "$BASHRC" || echo "$LINE" >> "$BASHRC"
done

export LD_LIBRARY_PATH="/usr/local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PV_PLUGIN_PATH="/usr/local/bin/plugins${PV_PLUGIN_PATH:+:$PV_PLUGIN_PATH}"
export PYTHONPATH="/usr/local/lib/python3.12/site-packages${PYTHONPATH:+:$PYTHONPATH}"
hash -r 2>/dev/null || true

echo
echo "Installation complete."
echo "Run ParaView with: paraview"
