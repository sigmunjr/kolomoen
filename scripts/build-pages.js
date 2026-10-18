"use strict";

// A static deployment: no Electron process, API server, or Cesium ion account.
const fs = require("node:fs");
const path = require("node:path");
const { copyCesium } = require("../desktop/copy-cesium");

const projectRoot = path.join(__dirname, "..");
const results = path.join(projectRoot, "results");
const output = path.join(projectRoot, "dist", "pages");
const data = JSON.parse(fs.readFileSync(path.join(results, "viewer_data.json"), "utf8"));

// Keep the public site limited to the viewer, its data, and the result downloads.
// Exclude the old 21 MB embedded HTML and the analysis-only GeoTIFF.
const files = [
  "index.html", "viewer_data.json", "summary.json", "critical_plateaus.geojson",
  "hillshade.png", data.terrain.raw_file,
  ...data.terrain.levels.map(level => level.file),
  "RESULTAT.md", "routes.geojson", "critical_plateaus.csv", "input_points.csv",
  "connectivity_checks.csv", "route_all_cells.csv", "route_50m.csv",
  "reference_27_NEW_points.csv", "reference_50m.csv", "reference_all_crossed_cells.csv",
  "valley_route_all_cells.csv", "saddle_north_south_2km.csv",
  "overview.png", "profiles.png", "saddle_section.png"
];

for (const file of files) {
  if (path.basename(file) !== file || !fs.statSync(path.join(results, file)).isFile()) {
    throw new Error(`Missing or invalid viewer asset: ${file}`);
  }
  const start = fs.readFileSync(path.join(results, file)).subarray(0, 100).toString();
  if (start.startsWith("version https://git-lfs.github.com/spec/v1")) {
    throw new Error(`Git LFS pointer instead of file contents: ${file}`);
  }
}
const expectedSize = data.terrain.width * data.terrain.height * 4;
if (fs.statSync(path.join(results, data.terrain.raw_file)).size !== expectedSize) {
  throw new Error(`terrain.f32 must contain ${expectedSize} bytes.`);
}

// Only remove this generated subdirectory, preserving desktop installer builds.
fs.rmSync(output, { recursive: true, force: true });
fs.mkdirSync(output, { recursive: true });
for (const file of files) fs.copyFileSync(path.join(results, file), path.join(output, file));
copyCesium(path.join(output, "vendor", "Cesium"));
fs.writeFileSync(path.join(output, ".nojekyll"), "");

// Catch broken local links before a workflow uploads an incomplete site.
const html = fs.readFileSync(path.join(output, "index.html"), "utf8");
for (const [, link] of html.matchAll(/(?:href|src)="([^"]+)"/g)) {
  if (/^(?:https?:|data:|#)/.test(link)) continue;
  if (link.startsWith("/") || !fs.existsSync(path.join(output, link))) {
    throw new Error(`Broken or root-relative site link: ${link}`);
  }
}
for (const directory of ["Workers", "Assets", "ThirdParty", "Widgets"]) {
  if (!fs.statSync(path.join(output, "vendor", "Cesium", directory)).isDirectory()) {
    throw new Error(`Missing Cesium ${directory} directory`);
  }
}
console.log(`GitHub Pages site ready: ${output}`);
