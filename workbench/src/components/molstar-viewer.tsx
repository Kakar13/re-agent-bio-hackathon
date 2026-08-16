"use client";

import { Boxes } from "lucide-react";
import { useMemo, useState } from "react";

type MolstarViewerProps = {
  structureUrl?: string;
  structurePath?: string;
  chainId?: string;
  residueIds?: string[];
  mappingStatus?: "verified_exact_sequence" | "verified_terminal_trim";
  unresolvedSequencePositions?: number[];
  spatialTracks?: Record<string, number[]>;
};

function scriptValue(value: unknown) {
  return JSON.stringify(value)
    .replaceAll("<", "\\u003c")
    .replaceAll("\u2028", "\\u2028")
    .replaceAll("\u2029", "\\u2029");
}

function preferredTrack(names: string[]) {
  return (
    names.find((name) => name === "mhci_processing_risk_max") ??
    names.find((name) => name === "mhci_presentation_propensity_mean") ??
    names.find((name) => name === "response_propensity") ??
    names.find((name) => name === "netmhciipan_el_support") ??
    names[0] ??
    ""
  );
}

export function MolstarViewer({
  structureUrl,
  structurePath,
  chainId = "A",
  residueIds = [],
  mappingStatus = "verified_exact_sequence",
  unresolvedSequencePositions = [],
  spatialTracks = {},
}: MolstarViewerProps) {
  const validTracks = useMemo(
    () =>
      Object.entries(spatialTracks).filter(
        ([, values]) =>
          values.length === residueIds.length &&
          values.length > 0 &&
          values.every((value) => Number.isFinite(value)),
      ),
    [residueIds.length, spatialTracks],
  );
  const trackNames = useMemo(() => validTracks.map(([name]) => name), [validTracks]);
  const defaultTrack = preferredTrack(trackNames);
  const [requestedTrack, setRequestedTrack] = useState("");
  const selectedTrack = trackNames.includes(requestedTrack) ? requestedTrack : defaultTrack;

  if (!structureUrl) {
    return (
      <div className="molstar-empty">
        <Boxes size={26} strokeWidth={1.4} />
        <strong>No validated structure attached</strong>
        <span>A sequence-matched AlphaFold2 PDB is required for the 3D evidence heatmap.</span>
      </div>
    );
  }

  const selectedValues =
    validTracks.find(([name]) => name === selectedTrack)?.[1] ?? validTracks[0]?.[1] ?? [];
  const hasHeatmap = selectedValues.length > 0;
  const peak = selectedValues.reduce((maximum, value) => Math.max(maximum, value), 0);
  const url = scriptValue(structureUrl);
  const chain = scriptValue(chainId);
  const ids = scriptValue(residueIds);
  const values = scriptValue(selectedValues);
  const document = `<!doctype html>
<html><head>
<meta charset="utf-8">
<link rel="stylesheet" href="/molstar/molstar.css">
<style>
html,body,#app{width:100%;height:100%;margin:0;overflow:hidden}
#status{position:absolute;inset:12px auto auto 12px;z-index:20;max-width:70%;padding:7px 9px;
font:12px/1.35 system-ui,sans-serif;color:#7f1d1d;background:#fef2f2;border:1px solid #fecaca;
display:none}
</style>
</head><body><div id="status"></div><div id="app"></div>
<script src="/molstar/molstar.js"></script>
<script>
function rewritePdbBFactors(pdb, targetChain, targetResidueIds, trackValues) {
  var residueIndex = new Map(targetResidueIds.map(function(id, index) { return [id, index]; }));
  return pdb.split(/\\r?\\n/).map(function(line) {
    if (!(line.startsWith("ATOM  ") || line.startsWith("HETATM"))) return line;
    var padded = line.padEnd(80, " ");
    if (padded.slice(21, 22).trim() !== targetChain) return null;
    var residueId = padded.slice(22, 26).trim() + padded.slice(26, 27).trim();
    var index = residueIndex.get(residueId);
    var value = index === undefined ? 0 : Math.max(0, Math.min(1, trackValues[index] || 0));
    // Mol* maps low B-factors to blue and high B-factors to red. Store the
    // selected [0,1] evidence value directly so its direction matches the UI.
    var riskScore = (100 * value).toFixed(2).padStart(6, " ");
    return padded.slice(0, 60) + riskScore + padded.slice(66);
  }).filter(function(line) { return line !== null; }).join("\\n");
}

async function start() {
  var response = await fetch(${url});
  if (!response.ok) throw new Error("Structure request failed with status " + response.status);
  var pdb = await response.text();
  var viewer = await molstar.Viewer.create("app", {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowSequence: true,
    layoutShowLog: false,
    viewportShowExpand: false
  });
  var trackValues = ${values};
  // Render only the scored chain, even when no residue-level track is available.
  pdb = rewritePdbBFactors(pdb, ${chain}, ${ids}, trackValues);
  await viewer.loadStructureFromData(pdb, "pdb");
  if (trackValues.length) {
    var structures = viewer.plugin.managers.structure.hierarchy.current.structures;
    var components = structures.flatMap(function(structure) { return structure.components; });
    await viewer.plugin.managers.structure.component.updateRepresentationsTheme(
      components,
      {
        color: "uncertainty",
        colorParams: {
          domain: [0, 100],
          list: {
            kind: "interpolate",
            colors: [0xdc2626, 0xf7f7f7, 0x1e40af]
          }
        }
      }
    );
  }
}

start().catch(function(error) {
  var status = document.getElementById("status");
  status.textContent = "Mol* could not render this structure: " + error.message;
  status.style.display = "block";
});
</script></body></html>`;

  return (
    <div className="molstar-shell">
      <div className="molstar-toolbar">
        <div>
          <span className="eyebrow">
            Showing scored chain {chainId} only
            {mappingStatus === "verified_terminal_trim"
              ? ` · ${residueIds.length - unresolvedSequencePositions.length}/${
                  residueIds.length
                } residues resolved`
              : ""}
          </span>
          <strong>{structurePath ?? "Attached PDB structure"}</strong>
        </div>
        {validTracks.length > 0 && (
          <label>
            3D evidence track
            <select
              aria-label="3D evidence track"
              value={selectedTrack || validTracks[0][0]}
              onChange={(event) => setRequestedTrack(event.target.value)}
            >
              {validTracks.map(([name]) => (
                <option key={name} value={name}>
                  {name.replaceAll("_", " ")}
                </option>
              ))}
            </select>
          </label>
        )}
      </div>
      <iframe
        key={`${structureUrl}-${selectedTrack}`}
        className="molstar-frame"
        title="Molstar 3D residue evidence heatmap"
        srcDoc={document}
        sandbox="allow-scripts allow-same-origin"
      />
      {hasHeatmap ? (
        <div className="molstar-legend">
          <span><i className="legend-low" />Low track value</span>
          <span>{selectedTrack.replaceAll("_", " ")} peak {peak.toFixed(3)}</span>
          <span><i className="legend-high" />High track value</span>
        </div>
      ) : (
        <p className="empty-copy">The structure is available, but no full-length residue track maps to it.</p>
      )}
    </div>
  );
}
