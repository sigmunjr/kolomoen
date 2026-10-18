"use strict";

const { app, BrowserWindow, shell } = require("electron");
const http = require("node:http");
const fs = require("node:fs");
const path = require("node:path");

const MIME_TYPES = {
  ".css": "text/css; charset=utf-8",
  ".csv": "text/csv; charset=utf-8",
  ".f32": "application/octet-stream",
  ".gif": "image/gif",
  ".geojson": "application/geo+json; charset=utf-8",
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".md": "text/markdown; charset=utf-8",
  ".png": "image/png",
  ".svg": "image/svg+xml",
  ".wasm": "application/wasm",
  ".webp": "image/webp",
  ".xml": "application/xml; charset=utf-8"
};

let server;

function viewerDirectory() {
  return app.isPackaged
    ? path.join(process.resourcesPath, "results")
    : path.join(__dirname, "..", "results");
}

function startLocalServer(root) {
  return new Promise((resolve, reject) => {
    server = http.createServer((request, response) => {
      try {
        const requestUrl = new URL(request.url, "http://127.0.0.1");
        let relativePath = decodeURIComponent(requestUrl.pathname).replace(/^\/+/, "");
        if (!relativePath) relativePath = "index.html";
        const filePath = path.resolve(root, relativePath);
        const relative = path.relative(root, filePath);
        if (relative.startsWith("..") || path.isAbsolute(relative)) {
          response.writeHead(403).end("Forbidden");
          return;
        }
        const stat = fs.statSync(filePath);
        if (!stat.isFile()) throw new Error("Not a file");
        response.writeHead(200, {
          "Content-Type": MIME_TYPES[path.extname(filePath).toLowerCase()] || "application/octet-stream",
          "Cache-Control": app.isPackaged ? "public, max-age=86400" : "no-cache",
          "Cross-Origin-Opener-Policy": "same-origin"
        });
        fs.createReadStream(filePath).pipe(response);
      } catch {
        response.writeHead(404, { "Content-Type": "text/plain; charset=utf-8" }).end("Not found");
      }
    });
    server.once("error", reject);
    server.listen(0, "127.0.0.1", () => resolve(server.address().port));
  });
}

async function createWindow() {
  const port = await startLocalServer(viewerDirectory());
  const window = new BrowserWindow({
    width: 1440,
    height: 920,
    minWidth: 900,
    minHeight: 650,
    backgroundColor: "#111b27",
    title: "Heradsbygd–Kolomoen",
    webPreferences: {
      contextIsolation: true,
      nodeIntegration: false,
      sandbox: true
    }
  });
  window.removeMenu();
  window.webContents.setWindowOpenHandler(({ url }) => {
    if (url.startsWith("http://") || url.startsWith("https://")) shell.openExternal(url);
    return { action: "deny" };
  });
  await window.loadURL(`http://127.0.0.1:${port}/`);
}

app.whenReady().then(createWindow).catch((error) => {
  console.error(error);
  app.quit();
});

app.on("window-all-closed", () => {
  if (server) server.close();
  app.quit();
});
