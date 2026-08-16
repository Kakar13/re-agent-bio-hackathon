"use client";

import clsx from "clsx";

export function HeatmapTrack({
  label,
  values,
  active = false,
}: {
  label: string;
  values: number[];
  active?: boolean;
}) {
  const peak = values.length ? Math.max(...values) : 0;

  return (
    <div className="track-row">
      <div className="track-label">
        <span className={clsx("track-dot", active && "track-dot-active")} />
        <span>{label.replaceAll("_", " ")}</span>
        <span className="track-peak">{peak.toFixed(2)}</span>
      </div>
      <div className="heatmap" aria-label={`${label} residue heatmap`}>
        {values.map((value, index) => (
          <span
            key={`${label}-${index}`}
            className="heat-cell"
            style={{
              opacity: 0.16 + Math.min(1, value) * 0.84,
              height: `${Math.max(18, value * 100)}%`,
            }}
            title={`Residue ${index + 1}: ${value.toFixed(3)}`}
          />
        ))}
      </div>
    </div>
  );
}
