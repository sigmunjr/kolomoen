"use strict";

const fs = require("node:fs");
const path = require("node:path");

const projectRoot = path.join(__dirname, "..");
const source = path.join(projectRoot, "node_modules", "cesium", "Build", "Cesium");
function copyCesium(destination = path.join(projectRoot, "results", "vendor", "Cesium")) {
  if (!fs.existsSync(source)) {
    throw new Error("Cesium is not installed. Run npm ci first.");
  }
  fs.mkdirSync(path.dirname(destination), { recursive: true });
  fs.cpSync(source, destination, { recursive: true });
  // Cesium's npm Build directory does not include the license at its root.
  fs.copyFileSync(path.join(projectRoot, "node_modules", "cesium", "LICENSE.md"),
    path.join(destination, "LICENSE.md"));
  console.log(`Copied local Cesium assets to ${destination}`);
}

if (require.main === module) copyCesium();

module.exports = { copyCesium };

