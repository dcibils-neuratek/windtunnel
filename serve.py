#!/usr/bin/env python3
"""Static dev server for the wind tunnel.

Two things `python -m http.server` will not do for us:

1. Never cache. Browsers hold on to .js and .css aggressively, which makes
   editing maddening — you reload and get the old file.

2. Cross-origin isolation (COOP + COEP). SharedArrayBuffer is gated behind
   it, and without SharedArrayBuffer the multi-core solver cannot run: the
   worker threads would have no way to share the lattice. This is also why
   multi-core is unavailable when index.html is opened straight off disk —
   a file:// URL has no headers at all.

    python serve.py [port]        # default 8777
"""
import sys
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class NoCacheHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        # unlocks SharedArrayBuffer -> the multi-core solver
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        super().end_headers()

    def log_message(self, fmt, *args):
        sys.stderr.write("%s\n" % (fmt % args))


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8777
    print("wind tunnel: http://localhost:%d  (no-cache)" % port)
    ThreadingHTTPServer(("127.0.0.1", port), NoCacheHandler).serve_forever()
