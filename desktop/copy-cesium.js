"use strict";

const fs = require("node:fs");
const path = require("node:path");

const projectRoot = path.join(__dirname, "..");
const source = path.join(projectRoot, "node_modules", "cesium", "Build", "Cesium");
const destination = path.join(projectRoot, "results", "vendor", "Cesium");

if (!fs.existsSync(source)) {
  throw new Error("Cesium is not installed. Run npm install first.");
}

fs.rmSync(destination, { recursive: true, force: true });
fs.mkdirSync(path.dirname(destination), { recursive: true });
fs.cpSync(source, destination, { recursive: true });
console.log(`Copied local Cesium assets to ${destination}`);

