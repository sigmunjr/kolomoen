#!/usr/bin/env python3
"""Serve the bundled CesiumJS results viewer on localhost and open it.

The viewer includes selectable ArcGIS World Imagery and OpenStreetMap background
layers. No third-party Python packages are needed to view existing results.
"""
from __future__ import annotations
import argparse
import functools
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import threading
import webbrowser


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error("Port must be between 1 and 65535.")
    directory = args.directory.resolve()
    if not (directory / "index.html").is_file():
        parser.error(f"No viewer found in {directory}. Run lowest_passage.py first.")
    handler = functools.partial(SimpleHTTPRequestHandler, directory=str(directory))
    try:
        server = ThreadingHTTPServer(("127.0.0.1", args.port), handler)
    except OSError as exc:
        parser.error(f"Cannot bind localhost:{args.port}: {exc}. Try --port 8001.")
    url = f"http://localhost:{args.port}/"
    print(f"Serving {directory}\nOpen {url}\nCtrl+C to stop.", flush=True)
    if not args.no_browser:
        threading.Timer(0.7, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
