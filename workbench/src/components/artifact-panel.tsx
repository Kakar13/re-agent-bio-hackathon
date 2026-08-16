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
import remarkGfm from "remark-gfm";

import { HeatmapTrack } from "@/components/heatmap-track";
import { MolstarViewer } from "@/components/molstar-viewer";
import type { Review, ScientificArtifact } from "@/lib/types";

type Tab =
  | "Overview"
  | "Objective"
  | "Literature"
  | "Toolchain"
  | "Candidates"
  | "Sequence"
  | "Structure"
  | "Immunogenicity"
  | "Evidence"
  | "Provenance"
  | "Review";

export function ArtifactPanel({
  artifact,
  review,
  onFork,
}: {
  artifact?: ScientificArtifact;
  review?: Review;
  onFork: () => void;
}) {
  const [tab, setTab] = useState<Tab>("Structure");
  const [candidateIndex, setCandidateIndex] = useState(0);
  const [annotations, setAnnotations] = useState<string[]>([]);
  const [draft, setDraft] = useState("");
  const architecture = artifact?.payload.architecture;
  const campaign = artifact?.payload.manifest;
  const assessments = artifact?.payload.assessments ?? [];
  const assessment = assessments[candidateIndex] ?? artifact?.payload.assessment;
  const structure = assessment?.structure;
  const structureUrl = structure
    ? `/api/structure?path=${encodeURIComponent(structure.path)}`
    : undefined;
  // Peak MHC-I processing risk is the one headline number this lane can stand
  // behind; the fused combined rank stays withheld without a response model.
  const mhciSummary = assessment?.mhc_i_surrogate_results.find((result) =>
    result.adapter_id.startsWith("team-e2e-pls-chao1"),
  )?.protein_summary;
  const peakMhciRisk = mhciSummary?.max_risk;
  const topMeanMhciRisk = mhciSummary?.top_k_mean_risk;

  const mhc = useMemo(
    () => assessment?.mhc_results.find((result) => result.provider_id === "netmhciipan"),
    [assessment],
  );
  const mhciSurrogates = assessment?.mhc_i_surrogate_results ?? [];
  const tabs = useMemo<Tab[]>(() => {
    if (!artifact) return ["Overview"];
    if (artifact.kind === "target_research" || artifact.kind === "paperclip_evidence") {
      return ["Overview", "Literature", "Provenance", "Review"];
    }
    if (artifact.kind === "proto_tool_validation") {
      return ["Overview", "Toolchain", "Provenance", "Review"];
    }
    if (artifact.kind === "design_spec" || artifact.kind === "design_campaign_plan") {
      return ["Overview", "Objective", "Toolchain", "Provenance", "Review"];
    }
    if (campaign) {
      return [
        "Overview",
        "Candidates",
        "Sequence",
        "Structure",
        "Immunogenicity",
        "Evidence",
        "Provenance",
        "Review",
      ];
    }
    return ["Overview", "Sequence", "Structure", "Evidence", "Provenance", "Review"];
  }, [artifact, campaign]);

  if (!artifact) {
    return (
      <aside className="artifact-panel">
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
        <div className="artifact-header-right">
          {peakMhciRisk != null && (
            <div className="headline-score" title="Peak MHC-I processing risk across all 9-mer windows">
              <span className="headline-score-value">{(peakMhciRisk * 100).toFixed(1)}</span>
              <span className="headline-score-label">
                MHC-I risk
                {topMeanMhciRisk != null && ` · top-5 ${(topMeanMhciRisk * 100).toFixed(1)}`}
              </span>
            </div>
          )}
          <button className="icon-button labeled" onClick={onFork} title="Fork this session">
            <GitFork size={15} />
            Fork
          </button>
        </div>
      </header>

      <nav className="artifact-tabs" aria-label="Artifact views">
        {tabs.map((item) => (
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
              {!architecture && !campaign && peakMhciRisk != null && (
                <Metric
                  label="MHC-I risk · peak 9-mer"
                  value={(peakMhciRisk * 100).toFixed(1)}
                />
              )}
              {!architecture && !campaign && topMeanMhciRisk != null && (
                <Metric
                  label="MHC-I risk · top-5 mean"
                  value={(topMeanMhciRisk * 100).toFixed(1)}
                />
              )}
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

        {tab === "Objective" && (
          <div className="stack">
            <section className="canvas-card">
              <span className="eyebrow">Natural-language objective</span>
              <h3>{artifact.payload.objective ?? artifact.title}</h3>
              <p className="empty-copy">
                {artifact.payload.claim_boundary ?? "No claim boundary recorded."}
              </p>
            </section>
            <section className="canvas-card code-card">
              <h3>Validated design specification</h3>
              <pre>{JSON.stringify(artifact.payload.spec ?? { spec_path: artifact.payload.spec_path }, null, 2)}</pre>
            </section>
          </div>
        )}

        {tab === "Literature" && (
          <div className="stack">
            <section className="canvas-card">
              <span className="eyebrow">Paperclip evidence</span>
              <h3>Line-pinned research record</h3>
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
              <h3>Retrieved evidence</h3>
              <pre>{JSON.stringify(artifact.payload.research ?? artifact.payload.result ?? {}, null, 2)}</pre>
            </section>
          </div>
        )}

        {tab === "Toolchain" && (
          <div className="stack">
            <section className="canvas-card">
              <span className="eyebrow">Application-owned execution</span>
              <h3>Proto tool contracts and phases</h3>
              <div className="lane-list">
                {(campaign?.executions ?? []).map((execution) => (
                  <Lane
                    key={execution.tool_key}
                    name={execution.tool_key}
                    status={execution.status}
                  />
                ))}
                {!campaign?.executions?.length &&
                  Object.keys(artifact.payload.schemas ?? {}).map((toolKey) => (
                    <Lane key={toolKey} name={toolKey} status="ready" />
                  ))}
              </div>
            </section>
            <section className="canvas-card code-card">
              <h3>Runtime record</h3>
              <pre>
                {JSON.stringify(
                  {
                    runtime: artifact.payload.runtime,
                    workspace: artifact.payload.workspace,
                    phases: campaign?.phase_status,
                  },
                  null,
                  2,
                )}
              </pre>
            </section>
          </div>
        )}

        {tab === "Candidates" && (
          <section className="canvas-card">
            <div className="section-heading">
              <div>
                <span className="eyebrow">Generated designs</span>
                <h3>Candidate comparison</h3>
              </div>
              <span className="status-chip">{campaign?.status ?? "missing"}</span>
            </div>
            <div className="data-table-wrap">
              <table className="data-table">
                <thead>
                  <tr>
                    <th>Candidate</th>
                    <th>Validation</th>
                    <th>Screening</th>
                    <th>Passed gates</th>
                  </tr>
                </thead>
                <tbody>
                  {(campaign?.candidates ?? []).map((candidate) => {
                    const assessmentIndex = assessments.findIndex(
                      (row) => row.candidate_id === candidate.candidate_id,
                    );
                    return (
                      <tr
                        key={candidate.candidate_id}
                        className={assessmentIndex === candidateIndex ? "selected-row" : undefined}
                        onClick={() => assessmentIndex >= 0 && setCandidateIndex(assessmentIndex)}
                      >
                        <td className="mono">{candidate.candidate_id}</td>
                        <td>{candidate.validation_status}</td>
                        <td>{candidate.screening_status}</td>
                        <td>
                          {candidate.validation_checks.filter((check) => check.passed).length}/
                          {candidate.validation_checks.length}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
            <p className="empty-copy">
              Select a screened candidate, then open Structure, Immunogenicity, or Evidence.
            </p>
          </section>
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

        {tab === "Immunogenicity" && (
          <div className="stack">
            <section className="canvas-card">
              <div className="section-heading">
                <div>
                  <span className="eyebrow">Selected candidate</span>
                  <h3>{assessment?.candidate_id ?? "No screened candidate"}</h3>
                </div>
                <ShieldCheck size={20} />
              </div>
              <div className="metric-grid">
                <Metric
                  label="Combined rank"
                  value={
                    assessment?.combined_rank_score == null
                      ? "Withheld"
                      : assessment.combined_rank_score.toFixed(3)
                  }
                  warning={assessment?.combined_rank_score == null}
                />
                <Metric label="MHC-II alleles" value={String(mhc?.supported_alleles.length ?? 0)} />
                <Metric label="EL + BA hits" value={String(mhc?.hits.length ?? 0)} />
                <Metric label="MHC-I lanes" value={String(mhciSurrogates.length)} />
              </div>
            </section>
            <section className="canvas-card">
              <h3>Independent evidence lanes</h3>
              <div className="lane-list">
                <Lane name="Response propensity" status={assessment?.response_results[0]?.status} />
                <Lane name="NetMHCIIpan EL" status={mhc?.status} />
                <Lane name="NetMHCIIpan BA" status={mhc?.status} />
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
          </div>
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
  const ok =
    status === "ok" ||
    status === "ready" ||
    status === "separate" ||
    status === "completed";
  return (
    <div className="lane">
      <span className={ok ? "lane-status ok" : "lane-status"} />
      <strong>{name}</strong>
      <span>{status ?? "missing"}</span>
    </div>
  );
}
