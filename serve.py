#!/usr/bin/env python3
"""Static dev server for the wind tunnel.

Identical to `python -m http.server` except it tells the browser never to
cache anything. Browsers are aggressive about holding on to .js and .css,
which makes editing this app maddening — you reload and get the old file.

    python serve.py [port]        # default 8777
"""
import sys
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class NoCacheHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

    def log_message(self, fmt, *args):
        sys.stderr.write("%s\n" % (fmt % args))


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8777
    print("wind tunnel: http://localhost:%d  (no-cache)" % port)
    ThreadingHTTPServer(("127.0.0.1", port), NoCacheHandler).serve_forever()
