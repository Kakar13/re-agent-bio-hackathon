"use client";

import {
  Bot,
  Check,
  CircleAlert,
  ClipboardCheck,
  Dna,
  FileCode2,
  GitFork,
  LoaderCircle,
  MessageSquarePlus,
  ShieldCheck,
} from "lucide-react";
import { useMemo, useState } from "react";
import ReactMarkdown from "react-markdown";

import { HeatmapTrack } from "@/components/heatmap-track";
import { MolstarViewer } from "@/components/molstar-viewer";
import type { Review, ScientificArtifact } from "@/lib/types";

const TABS = ["Overview", "Sequence", "Structure", "Evidence", "Provenance", "Review"] as const;
type Tab = (typeof TABS)[number];

export function ArtifactPanel({
  artifact,
  review,
  onFork,
  agentPrompt,
  agentResponse,
  agentIsLoading,
}: {
  artifact?: ScientificArtifact;
  review?: Review;
  onFork: () => void;
  agentPrompt: string;
  agentResponse: string;
  agentIsLoading: boolean;
}) {
  const [tab, setTab] = useState<Tab>("Overview");
  const [annotations, setAnnotations] = useState<string[]>([]);
  const [draft, setDraft] = useState("");
  const assessment = artifact?.payload.assessment;
  const architecture = artifact?.payload.architecture;
  const campaign = artifact?.payload.manifest;
  const structure = assessment?.structure;
  const structureUrl = structure
    ? `/api/structure?path=${encodeURIComponent(structure.path)}`
    : undefined;

  const mhc = useMemo(
    () => assessment?.mhc_results.find((result) => result.provider_id === "netmhciipan"),
    [assessment],
  );
  const mhciSurrogates = assessment?.mhc_i_surrogate_results ?? [];

  if (!artifact) {
    return (
      <aside className="artifact-panel">
        <AgentResponseDock
          prompt={agentPrompt}
          response={agentResponse}
          isLoading={agentIsLoading}
        />
        <div className="artifact-empty">
          <FileCode2 size={30} strokeWidth={1.4} />
          <h2>No artifact selected</h2>
          <p>Run an architecture check or candidate screen to populate this canvas.</p>
        </div>
      </aside>
    );
  }

  function addAnnotation() {
    const value = draft.trim();
    if (!value) return;
    setAnnotations((current) => [...current, value]);
    setDraft("");
  }

  return (
    <aside className="artifact-panel">
      <header className="artifact-header">
        <div>
          <div className="eyebrow">{artifact.kind.replaceAll("_", " ")}</div>
          <h2>{artifact.title}</h2>
          <p className="artifact-path">{artifact.path}</p>
        </div>
        <button className="icon-button labeled" onClick={onFork} title="Fork this session">
          <GitFork size={15} />
          Fork
        </button>
      </header>

      <AgentResponseDock
        prompt={agentPrompt}
        response={agentResponse}
        isLoading={agentIsLoading}
      />

      <nav className="artifact-tabs" aria-label="Artifact views">
        {TABS.map((item) => (
          <button
            key={item}
            className={item === tab ? "active" : undefined}
            onClick={() => setTab(item)}
          >
            {item}
          </button>
        ))}
      </nav>

      <div className="artifact-content">
        {tab === "Overview" && (
          <div className="stack">
            <section className="claim-boundary">
              <ShieldCheck size={18} />
              <div>
                <span>Claim boundary</span>
                <p>{artifact.payload.claim_boundary ?? "No claim boundary recorded."}</p>
              </div>
            </section>
            <div className="metric-grid">
              <Metric
                label={architecture ? "Architecture" : campaign ? "Candidates" : "Combined rank"}
                value={
                  architecture
                    ? architecture.status
                    : campaign
                    ? `${campaign.candidate_counts.screened ?? 0} / ${
                        campaign.candidate_counts.generated ?? 0
                      } screened`
                    : assessment?.combined_rank_score == null
                    ? "Withheld"
                    : assessment.combined_rank_score.toFixed(3)
                }
                warning={!architecture && assessment?.combined_rank_score == null}
              />
              <Metric
                label={campaign ? "Validation passed" : "MHC-II alleles"}
                value={
                  campaign
                    ? String(campaign.candidate_counts.validation_passed ?? 0)
                    : architecture
                    ? "18 planned"
                    : String(mhc?.supported_alleles.length ?? 0)
                }
              />
              <Metric
                label={architecture ? "MHC channels" : "EL + BA hits"}
                value={
                  architecture
                    ? architecture.providers.netmhciipan.channels.join(" + ")
                    : String(mhc?.hits.length ?? 0)
                }
              />
              <Metric
                label="Reviewer"
                value={review?.status ?? "Pending"}
                warning={review?.status === "fail"}
              />
            </div>
            {campaign && (
              <section className="canvas-card">
                <div className="section-heading">
                  <div>
                    <span className="eyebrow">Structural validation</span>
                    <h3>Candidate handoff gates</h3>
                  </div>
                  <span className="status-chip">{campaign.status}</span>
                </div>
                <div className="data-table-wrap">
                  <table className="data-table">
                    <thead>
                      <tr>
                        <th>Candidate</th>
                        <th>Validation</th>
                        <th>Screening</th>
                        <th>Checks</th>
                      </tr>
                    </thead>
                    <tbody>
                      {campaign.candidates.slice(0, 12).map((candidate) => (
                        <tr key={candidate.candidate_id}>
                          <td className="mono">{candidate.candidate_id}</td>
                          <td>{candidate.validation_status}</td>
                          <td>{candidate.screening_status}</td>
                          <td>
                            {candidate.validation_checks.filter((check) => check.passed).length}/
                            {candidate.validation_checks.length}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </section>
            )}
            <section className="canvas-card">
              <h3>Evidence lanes</h3>
              <div className="lane-list">
                <Lane
                  name="Response propensity"
                  status={
                    architecture
                      ? architecture.providers.response_model.available
                        ? "ready"
                        : "unavailable"
                      : assessment?.response_results[0]?.status
                  }
                />
                <Lane name="NetMHCIIpan EL" status={architecture ? "ready" : mhc?.status} />
                <Lane name="NetMHCIIpan BA" status={architecture ? "ready" : mhc?.status} />
                {architecture && (
                  <Lane name="Custom MHC-I processing models" status="separate" />
                )}
                {mhciSurrogates.map((surrogate) => (
                  <Lane
                    key={surrogate.adapter_id}
                    name={`${customModelLabel(surrogate.adapter_id)} MHC-I surrogate`}
                    status={surrogate.status}
                  />
                ))}
                <Lane name="Processing" status={assessment?.processing ? "ok" : "unavailable"} />
                <Lane name="Tolerance" status={assessment?.tolerance ? "ok" : "unavailable"} />
              </div>
            </section>
            {!!assessment?.warnings.length && (
              <section className="warning-list">
                {assessment.warnings.map((warning) => (
                  <p key={warning}>
                    <CircleAlert size={15} />
                    {warning}
                  </p>
                ))}
              </section>
            )}
          </div>
        )}

        {tab === "Sequence" && (
          <section className="canvas-card sequence-card">
            <div className="section-heading">
              <div>
                <span className="eyebrow">Candidate sequence</span>
                <h3>{assessment?.candidate_id ?? artifact.id}</h3>
              </div>
              <Dna size={21} />
            </div>
            <div className="sequence-ruler">
              {(assessment?.sequence ?? "No sequence attached").split("").map((residue, index) => (
                <span key={`${residue}-${index}`} title={`Residue ${index + 1}`}>
                  <b>{residue}</b>
                  {(index + 1) % 10 === 0 && <small>{index + 1}</small>}
                </span>
              ))}
            </div>
          </section>
        )}

        {tab === "Structure" && (
          <MolstarViewer
            structureUrl={structureUrl}
            structurePath={structure?.path}
            chainId={structure?.chain_id}
            residueIds={structure?.residue_ids}
            mappingStatus={structure?.mapping_status}
            unresolvedSequencePositions={structure?.unresolved_sequence_positions}
            spatialTracks={assessment?.spatial_tracks}
          />
        )}

        {tab === "Evidence" && (
          <div className="stack">
            <section className="canvas-card">
              <div className="section-heading">
                <div>
                  <span className="eyebrow">Residue-space projection</span>
                  <h3>Spatial evidence tracks</h3>
                </div>
                <span className="mono-label">{assessment?.sequence.length ?? 0} aa</span>
              </div>
              <div className="tracks">
                {Object.entries(assessment?.spatial_tracks ?? {}).map(([label, values]) => (
                  <HeatmapTrack key={label} label={label} values={values} active={label.includes("el")} />
                ))}
                {!Object.keys(assessment?.spatial_tracks ?? {}).length && (
                  <p className="empty-copy">No spatial tracks recorded.</p>
                )}
              </div>
            </section>
            <section className="canvas-card">
              <div className="section-heading">
                <div>
                  <span className="eyebrow">Independent MHC channels</span>
                  <h3>Top presentation and binding rows</h3>
                </div>
                <span className="status-chip">{mhc?.status ?? "missing"}</span>
              </div>
              <div className="data-table-wrap">
                <table className="data-table">
                  <thead>
                    <tr>
                      <th>Allele</th>
                      <th>Core</th>
                      <th>EL rank</th>
                      <th>BA rank</th>
                      <th>IC50 nM</th>
                    </tr>
                  </thead>
                  <tbody>
                    {(mhc?.hits ?? []).slice(0, 12).map((hit, index) => (
                      <tr key={`${hit.allele}-${hit.start}-${index}`}>
                        <td>{hit.allele}</td>
                        <td className="mono">{hit.core}</td>
                        <td>{hit.el_rank?.toFixed(2) ?? "—"}</td>
                        <td>{hit.ba_rank?.toFixed(2) ?? "—"}</td>
                        <td>{hit.ba_ic50_nm?.toFixed(0) ?? "—"}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </section>
            {mhciSurrogates.map((surrogate) => (
              <section className="canvas-card" key={surrogate.adapter_id}>
                <div className="section-heading">
                  <div>
                    <span className="eyebrow">Independent MHC-I processing lane</span>
                    <h3>Top {customModelLabel(surrogate.adapter_id)} 9-mers</h3>
                  </div>
                  <span className="status-chip">{surrogate.status}</span>
                </div>
                <div className="data-table-wrap">
                  <table className="data-table">
                    <thead>
                      <tr>
                        <th>Peptide</th>
                        <th>Cleavage N</th>
                        <th>Cleavage C</th>
                        <th>MHC-I</th>
                        <th>Composite</th>
                      </tr>
                    </thead>
                    <tbody>
                      {[...surrogate.predictions]
                        .sort(
                          (left, right) =>
                            right.composite_processing_risk -
                            left.composite_processing_risk,
                        )
                        .slice(0, 12)
                        .map((prediction) => (
                          <tr key={`${prediction.start}-${prediction.peptide}`}>
                            <td className="mono">{prediction.peptide}</td>
                            <td>{prediction.cleavage_n_probability.toFixed(3)}</td>
                            <td>{prediction.cleavage_c_probability.toFixed(3)}</td>
                            <td>{prediction.mhc_i_presentation_propensity.toFixed(3)}</td>
                            <td>{prediction.composite_processing_risk.toFixed(3)}</td>
                          </tr>
                        ))}
                    </tbody>
                  </table>
                </div>
                <p className="empty-copy">
                  HLA-A*02:01 MHC-I surrogate only; excluded from MHC-II fusion.
                </p>
              </section>
            ))}
          </div>
        )}

        {tab === "Provenance" && (
          <div className="stack">
            <section className="canvas-card">
              <h3>Immutable artifact</h3>
              <dl className="provenance-list">
                <div><dt>Run ID</dt><dd>{artifact.id}</dd></div>
                <div><dt>SHA-256</dt><dd className="mono">{artifact.sha256}</dd></div>
                <div><dt>Manifest</dt><dd>{artifact.manifest_path}</dd></div>
              </dl>
            </section>
            <section className="canvas-card">
              <h3>Method citations</h3>
              <div className="citation-list">
                {(artifact.payload.citations ?? []).map((citation) => (
                  <a key={citation.url} href={citation.url} target="_blank" rel="noreferrer">
                    <FileCode2 size={15} />
                    <span>{citation.claim}</span>
                  </a>
                ))}
              </div>
            </section>
            <section className="canvas-card code-card">
              <h3>Provider record</h3>
              <pre>{JSON.stringify(mhc?.provenance ?? {}, null, 2)}</pre>
            </section>
          </div>
        )}

        {tab === "Review" && (
          <div className="stack">
            <section className="canvas-card">
              <div className="section-heading">
                <div>
                  <span className="eyebrow">Deterministic gate</span>
                  <h3>{review?.summary ?? "Reviewer has not run"}</h3>
                </div>
                <ClipboardCheck size={21} />
              </div>
              <div className="review-checks">
                {(review?.checks ?? []).map((check) => (
                  <div key={check.name} className={check.passed ? "pass" : "fail"}>
                    {check.passed ? <Check size={15} /> : <CircleAlert size={15} />}
                    <div>
                      <strong>{check.name.replaceAll("_", " ")}</strong>
                      <p>{check.detail}</p>
                    </div>
                  </div>
                ))}
              </div>
            </section>
            <section className="canvas-card annotation-card">
              <div className="section-heading">
                <div>
                  <span className="eyebrow">Human review</span>
                  <h3>Inline annotations</h3>
                </div>
                <MessageSquarePlus size={20} />
              </div>
              <div className="annotation-input">
                <input
                  value={draft}
                  onChange={(event) => setDraft(event.target.value)}
                  onKeyDown={(event) => event.key === "Enter" && addAnnotation()}
                  placeholder="Add a reviewer note…"
                />
                <button onClick={addAnnotation}>Add</button>
              </div>
              {annotations.map((annotation, index) => (
                <p className="annotation" key={`${annotation}-${index}`}>{annotation}</p>
              ))}
            </section>
          </div>
        )}
      </div>
    </aside>
  );
}

function AgentResponseDock({
  prompt,
  response,
  isLoading,
}: {
  prompt: string;
  response: string;
  isLoading: boolean;
}) {
  return (
    <section className="agent-response-dock" aria-label="Agent response">
      <header>
        <span><Bot size={14} /> Agent response</span>
        <span className={isLoading ? "agent-response-status running" : "agent-response-status"}>
          {isLoading && <LoaderCircle className="spin" size={12} />}
          {isLoading ? "Responding" : response ? "Complete" : "Ready"}
        </span>
      </header>
      {prompt && (
        <p className="agent-response-prompt" title={prompt}>
          <b>Responding to</b>
          {prompt}
        </p>
      )}
      <div className="agent-response-copy" aria-live="polite">
        {response ? (
          <ReactMarkdown>{response}</ReactMarkdown>
        ) : isLoading ? (
          <p>The evidence is already rendering. The agent interpretation will appear here.</p>
        ) : (
          <p>Run a screen or ask the agent a question to see its interpretation here.</p>
        )}
      </div>
    </section>
  );
}

function customModelLabel(adapterId: string) {
  return adapterId.replace("team-e2e-pls-", "");
}

function Metric({ label, value, warning = false }: { label: string; value: string; warning?: boolean }) {
  return (
    <div className="metric">
      <span>{label}</span>
      <strong className={warning ? "warning" : undefined}>{value}</strong>
    </div>
  );
}

function Lane({ name, status }: { name: string; status?: string }) {
  const ok = status === "ok" || status === "ready" || status === "separate";
  return (
    <div className="lane">
      <span className={ok ? "lane-status ok" : "lane-status"} />
      <strong>{name}</strong>
      <span>{status ?? "missing"}</span>
    </div>
  );
}
