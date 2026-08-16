"use client";

import {
  ArrowRight,
  CheckCircle2,
  Dna,
  Fingerprint,
  Flame,
  LoaderCircle,
  Play,
  ShieldCheck,
} from "lucide-react";
import { FormEvent, ReactNode, useState } from "react";

import { MolstarViewer } from "@/components/molstar-viewer";
import type { Assessment, ScientificArtifact } from "@/lib/types";

const PDA_EXAMPLE_SEQUENCE =
  "GHMDEALALAARAREVRPRALARYRELTDDEEEVAEVERMADLICAQRLPPEWVIQLLKEILEEVKANPEKADEMIEENRDDVMLRTLWVLATEIYRAP";

export function Chao1Screen({
  sequence,
  assessment,
  artifact,
  isLoading,
  onSequenceChange,
  onRun,
}: {
  sequence: string;
  assessment?: Assessment;
  artifact?: ScientificArtifact;
  isLoading: boolean;
  onSequenceChange: (sequence: string) => void;
  onRun: (sequence: string, candidateId: string) => void;
}) {
  const [error, setError] = useState("");
  const cleaned = sequence.replace(/\s+/g, "").toUpperCase();
  const surrogate = assessment?.mhc_i_surrogate_results.find(
    (result) => result.adapter_id === "team-e2e-pls-chao1",
  );
  const riskTrack = surrogate?.spatial_tracks.mhci_processing_risk_max ?? [];
  const summary = surrogate?.protein_summary;
  const maxRisk = summary?.max_risk ?? Math.max(0, ...riskTrack);
  const topMean = summary?.top_k_mean_risk ?? 0;
  const checkpointHash = surrogate?.provenance.parameters.checkpoint_sha256;
  const structure = assessment?.structure;
  const structureUrl = structure
    ? `/api/structure?path=${encodeURIComponent(structure.path)}`
    : undefined;
  const topWindows = [...(surrogate?.predictions ?? [])]
    .sort(
      (left, right) =>
        right.composite_processing_risk - left.composite_processing_risk,
    )
    .slice(0, 5);

  function submit(event: FormEvent) {
    event.preventDefault();
    if (!/^[ACDEFGHIKLMNPQRSTVWY]{15,}$/.test(cleaned)) {
      setError("Use at least 15 standard amino-acid letters.");
      return;
    }
    setError("");
    onRun(
      cleaned,
      cleaned === PDA_EXAMPLE_SEQUENCE
        ? "pda:9s14:0"
        : `visual-chao1-${Date.now().toString(36)}`,
    );
  }

  return (
    <section className="chao1-workspace">
      <header className="chao1-titlebar">
        <div>
          <span className="eyebrow">Primary workflow</span>
          <h2>Sequence risk map</h2>
        </div>
        <span className="model-verified">
          <ShieldCheck size={14} />
          Chao1
        </span>
      </header>

      <form className="sequence-screen-form" onSubmit={submit}>
        <div className="sequence-input-heading">
          <label htmlFor="chao1-sequence">
            <Dna size={15} />
            Protein sequence
          </label>
          <span>{cleaned.length} aa</span>
        </div>
        <textarea
          id="chao1-sequence"
          aria-label="Protein sequence"
          value={sequence}
          onChange={(event) => onSequenceChange(event.target.value.toUpperCase())}
          placeholder="Paste amino-acid sequence…"
          spellCheck={false}
          rows={3}
        />
        <div className="sequence-actions">
          <button
            type="button"
            className="example-button"
            onClick={() => {
              setError("");
              onSequenceChange(PDA_EXAMPLE_SEQUENCE);
            }}
          >
            Load PDA example
          </button>
          {error && <span className="sequence-error">{error}</span>}
          <button type="submit" className="run-chao1" disabled={isLoading}>
            {isLoading ? (
              <LoaderCircle className="spin" size={15} />
            ) : (
              <Play size={14} fill="currentColor" />
            )}
            {isLoading ? (surrogate ? "Agent responding…" : "Running model…") : "Run chao1"}
          </button>
        </div>
      </form>

      <div className="model-flow" aria-label="Chao1 screening flow">
        <FlowStep icon={<Dna size={15} />} label="Sequence" detail={`${cleaned.length || "—"} aa`} />
        <ArrowRight size={14} />
        <FlowStep
          icon={
            isLoading && !surrogate
              ? <LoaderCircle className="spin" size={15} />
              : <Fingerprint size={15} />
          }
          label="Chao1 model"
          detail={surrogate?.status === "ok" ? "Verified" : isLoading ? "Running" : "Ready"}
          active={isLoading && !surrogate}
        />
        <ArrowRight size={14} />
        <FlowStep
          icon={<Flame size={15} />}
          label="Risk map"
          detail={surrogate?.status === "ok" ? "Complete" : "Waiting"}
        />
      </div>

      {!assessment && !isLoading && (
        <div className="visual-empty">
          <div className="visual-empty-map" aria-hidden="true">
            {Array.from({ length: 48 }, (_, index) => (
              <span key={index} style={{ opacity: 0.14 + ((index * 17) % 70) / 100 }} />
            ))}
          </div>
          <strong>Paste a sequence or load the PDA example.</strong>
          <span>The result appears here as a residue-level risk map.</span>
        </div>
      )}

      {isLoading && !assessment && (
        <div className="visual-loading">
          <LoaderCircle className="spin" size={24} />
          <div>
            <strong>Running the actual chao1 checkpoint</strong>
            <span>ESM-2 encoding → cleavage/TAP/MHC-I heads → residue projection</span>
          </div>
        </div>
      )}

      {assessment && surrogate && (
        <div className="risk-result">
          <div className="risk-summary">
            <div
              className="risk-dial"
              style={{
                background: `conic-gradient(#b65432 ${maxRisk * 360}deg, #e9e5dc 0deg)`,
              }}
            >
              <div>
                <strong>{(maxRisk * 100).toFixed(1)}</strong>
                <span>/ 100</span>
              </div>
            </div>
            <div className="risk-summary-copy">
              <span className="eyebrow">Highest 9-mer score</span>
              <h3>Chao1 processing risk</h3>
              <div className="risk-mini-metrics">
                <span>
                  <b>{(topMean * 100).toFixed(1)}</b>
                  top-5 mean
                </span>
                <span>
                  <b>{summary?.n_windows ?? surrogate.predictions.length}</b>
                  windows
                </span>
                <span>
                  <b>{((summary?.mean_confidence ?? 0) * 100).toFixed(0)}%</b>
                  confidence
                </span>
              </div>
            </div>
            <div className="checkpoint-proof">
              <CheckCircle2 size={16} />
              <div>
                <strong>Actual checkpoint loaded</strong>
                <span>{surrogate.adapter_id}</span>
                <code title={checkpointHash}>{checkpointHash?.slice(0, 16)}…</code>
              </div>
            </div>
          </div>

          {structure && (
            <section className="main-structure-card">
              <div className="risk-section-heading">
                <div>
                  <span className="eyebrow">Experimental coordinates</span>
                  <h3>3D structure · PDB 9S14</h3>
                </div>
                <span className="structure-resolution">
                  {structure.residue_ids.length -
                    structure.unresolved_sequence_positions.length}
                  /{structure.residue_ids.length} residues resolved
                </span>
              </div>
              <MolstarViewer
                structureUrl={structureUrl}
                structurePath={structure.path}
                chainId={structure.chain_id}
                residueIds={structure.residue_ids}
                mappingStatus={structure.mapping_status}
                unresolvedSequencePositions={structure.unresolved_sequence_positions}
                spatialTracks={assessment.spatial_tracks}
              />
            </section>
          )}

          <section className="residue-map-card">
            <div className="risk-section-heading">
              <div>
                <span className="eyebrow">Residue projection</span>
                <h3>{assessment.candidate_id}</h3>
              </div>
              <div className="risk-legend">
                <span>Lower</span>
                <i />
                <span>Higher</span>
              </div>
            </div>
            <div className="residue-risk-map" aria-label="Chao1 residue risk heatmap">
              {assessment.sequence.split("").map((residue, index) => {
                const risk = riskTrack[index] ?? 0;
                return (
                  <span
                    key={`${residue}-${index}`}
                    className="risk-residue"
                    style={{
                      backgroundColor: riskColor(risk),
                      color: risk > 0.34 && risk < 0.66 ? "#26302d" : "#ffffff",
                    }}
                    title={`Residue ${index + 1} ${residue}: ${risk.toFixed(3)}`}
                  >
                    <b>{residue}</b>
                    <small>{index + 1}</small>
                  </span>
                );
              })}
            </div>
          </section>

          <section className="hotspot-card">
            <div className="risk-section-heading">
              <div>
                <span className="eyebrow">Highest scoring regions</span>
                <h3>Top chao1 9-mers</h3>
              </div>
              <span className="model-output-id">{artifact?.id}</span>
            </div>
            <div className="hotspot-bars">
              {topWindows.map((window) => (
                <div className="hotspot-row" key={`${window.start}-${window.peptide}`}>
                  <code>{window.peptide}</code>
                  <span>{window.start + 1}–{window.end}</span>
                  <i>
                    <b style={{ width: `${window.composite_processing_risk * 100}%` }} />
                  </i>
                  <strong>{(window.composite_processing_risk * 100).toFixed(1)}</strong>
                </div>
              ))}
            </div>
          </section>
        </div>
      )}
    </section>
  );
}

function riskColor(value: number) {
  const risk = Math.max(0, Math.min(1, value));
  const low = [30, 64, 175];
  const middle = [247, 247, 247];
  const high = [220, 38, 38];
  const start = risk <= 0.5 ? low : middle;
  const end = risk <= 0.5 ? middle : high;
  const amount = risk <= 0.5 ? risk * 2 : (risk - 0.5) * 2;
  const channels = start.map((channel, index) =>
    Math.round(channel + (end[index] - channel) * amount),
  );
  return `rgb(${channels.join(" ")})`;
}

function FlowStep({
  icon,
  label,
  detail,
  active = false,
}: {
  icon: ReactNode;
  label: string;
  detail: string;
  active?: boolean;
}) {
  return (
    <div className={`flow-step ${active ? "active" : ""}`}>
      {icon}
      <div>
        <strong>{label}</strong>
        <span>{detail}</span>
      </div>
    </div>
  );
}
