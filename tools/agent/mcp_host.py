"""A minimal MCP host for recording a chat with `kofl mcp`: it starts the
server once, keeps it open, and lets short commands send tool calls to it
through a Unix socket. Every call and answer goes into a transcript (JSON
lines), next to the journal the server writes itself.

    python3 tools/agent/mcp_host.py serve SOCKET TRANSCRIPT -- kofl mcp --read-root IN --write-root OUT --journal J
    python3 tools/agent/mcp_host.py call SOCKET TOOL '{"json": "arguments"}' [--approve yes|no] [--full]
    python3 tools/agent/mcp_host.py tools SOCKET
    python3 tools/agent/mcp_host.py stop SOCKET

The host declares the `elicitation` capability, so a build over the MCP
limits asks before it runs (AG-15). The person's answer is given with the
call: `--approve yes` accepts, `--approve no` declines; without it the host
cancels the question, and the build is refused (TOOL010). The question and
the answer are written into the transcript.
"""
import json
import os
import socket
import subprocess
import sys
import threading
import time


def serve(sock_path, transcript, argv):
    p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
    lock = threading.Lock()
    call_lock = threading.Lock()
    waiting = {}
    state = {"n": 0, "approve": None, "asked": None}
    log = open(transcript, "a", buffering=1)

    def write(m):
        with lock:
            p.stdin.write(json.dumps(m) + "\n")
            p.stdin.flush()

    def reader():
        for line in p.stdout:
            m = json.loads(line)
            if m.get("method") == "elicitation/create":
                state["asked"] = m["params"].get("message")
                a = state["approve"]
                if a is None:
                    res = {"action": "cancel"}
                elif a:
                    res = {"action": "accept", "content": {"approve": True}}
                else:
                    res = {"action": "decline"}
                state["answered"] = res["action"]
                write({"jsonrpc": "2.0", "id": m["id"], "result": res})
            elif "id" in m and "method" not in m and m["id"] in waiting:
                waiting[m["id"]]["answer"] = m
                waiting[m["id"]]["event"].set()

    threading.Thread(target=reader, daemon=True).start()

    def post(method, params):
        with lock:
            state["n"] += 1
            i = state["n"]
        w = {"event": threading.Event(), "answer": None}
        waiting[i] = w
        write({"jsonrpc": "2.0", "id": i, "method": method, "params": params})
        return w

    w = post("initialize", {"protocolVersion": "2025-11-25", "capabilities": {"elicitation": {}},
                            "clientInfo": {"name": "koflcraft-mcp-host", "version": "1"}})
    w["event"].wait()
    write({"jsonrpc": "2.0", "method": "notifications/initialized"})

    def handle(conn):
        req = json.loads(conn.makefile().readline())
        op = req["op"]
        if op == "call":
            with call_lock:
                state["approve"] = req.get("approve")
                state["asked"] = None
                t0 = time.time()
                w = post("tools/call", {"name": req["tool"], "arguments": req["arguments"]})
                w["event"].wait()
                ans = w["answer"]
                r = ans.get("result") or {}
                content = r.get("content") or [{"text": json.dumps(ans.get("error"))}]
                out = {"tool": req["tool"], "arguments": req["arguments"],
                       "seconds": round(time.time() - t0, 2),
                       "headline": content[0]["text"].split("\n")[0],
                       "envelope": r.get("structuredContent")}
                if state["asked"] is not None:
                    out["asked"] = state["asked"]
                    out["answered"] = state["answered"]
                log.write(json.dumps(out) + "\n")
            conn.sendall((json.dumps(out) + "\n").encode())
        elif op == "tools":
            w = post("tools/list", {})
            w["event"].wait()
            conn.sendall((json.dumps(w["answer"]["result"]) + "\n").encode())
        elif op == "stop":
            conn.sendall(b'{"ok": true}\n')
            p.stdin.close()
            p.wait()
            os._exit(0)
        conn.close()

    if os.path.exists(sock_path):
        os.remove(sock_path)
    s = socket.socket(socket.AF_UNIX)
    s.bind(sock_path)
    s.listen(8)
    while True:
        c, _ = s.accept()
        threading.Thread(target=handle, args=(c,), daemon=True).start()


def client(sock_path, req):
    s = socket.socket(socket.AF_UNIX)
    s.connect(sock_path)
    s.sendall((json.dumps(req) + "\n").encode())
    data = b""
    while not data.endswith(b"\n"):
        chunk = s.recv(1 << 16)
        if not chunk:
            break
        data += chunk
    return json.loads(data)


def show(out, full=False):
    env = out["envelope"] or {}
    print(f'{out["tool"]}: {env.get("status")} (exit {env.get("exit_code")}, {out["seconds"]} s)')
    if "asked" in out:
        print(f"  asked: {out['asked']}\n  answered: {out['answered']}")
    res = json.dumps(env.get("result"), indent=1 if full else None)
    print(res if full or len(res) < 3000 else res[:3000] + " ...")
    for d in env.get("diagnostics", []):
        print(f"  {d.get('code')}: {d.get('message')}")
    for lim in env.get("limitations", []):
        print(f"  limitation: {lim}")


def main():
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    if a[0] == "serve":
        i = a.index("--")
        serve(a[1], a[2], a[i + 1:])
    elif a[0] == "call":
        approve = None
        if "--approve" in a:
            approve = a[a.index("--approve") + 1] == "yes"
        args = json.loads(a[3]) if len(a) > 3 and not a[3].startswith("--") else {}
        show(client(a[1], {"op": "call", "tool": a[2], "arguments": args, "approve": approve}), "--full" in a)
    elif a[0] == "tools":
        print(" ".join(t["name"] for t in client(a[1], {"op": "tools"})["tools"]))
    elif a[0] == "stop":
        client(a[1], {"op": "stop"})


if __name__ == "__main__":
    main()
