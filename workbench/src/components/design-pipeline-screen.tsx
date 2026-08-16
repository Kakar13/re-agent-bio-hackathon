"use client";

import {
  BookOpen,
  Check,
  ChevronRight,
  Cpu,
  Dna,
  FlaskConical,
  LoaderCircle,
  Play,
  RotateCcw,
  ShieldCheck,
  Sparkles,
  Square,
} from "lucide-react";
import { FormEvent, useMemo, useState } from "react";

import type { ScientificArtifact } from "@/lib/types";

const EXAMPLE_OBJECTIVE =
  "Design an 80-residue protein binder to IL-7Rα at the cytokine-binding interface. " +
  "Avoid cysteine, preserve high predicted monomer confidence, and prioritize low " +
  "predicted immunogenicity after structural validation.";

const PIPELINE_KINDS = new Set([
  "target_research",
  "paperclip_evidence",
  "proto_tool_validation",
  "design_spec",
  "design_campaign_plan",
  "design_to_screen_campaign",
  "design_to_screen_campaign_replay",
]);

export type PipelineToolEvent = {
  id: string;
  name: string;
  completed: boolean;
  error?: boolean;
};

export function DesignPipelineScreen({
  objective,
  artifacts,
  toolEvents,
  isLoading,
  onObjectiveChange,
  onRun,
  onReplay,
  onStop,
  onSelectArtifact,
}: {
  objective: string;
  artifacts: ScientificArtifact[];
  toolEvents: PipelineToolEvent[];
  isLoading: boolean;
  onObjectiveChange: (value: string) => void;
  onRun: (objective: string) => void;
  onReplay: () => void;
  onStop: () => void;
  onSelectArtifact: (id: string) => void;
}) {
  const [error, setError] = useState("");
  const pipelineArtifacts = artifacts.filter((artifact) => PIPELINE_KINDS.has(artifact.kind));
  const latestCampaign = [...pipelineArtifacts]
    .reverse()
    .find((artifact) => artifact.payload.manifest);
  const phases = latestCampaign?.payload.manifest?.phase_status ?? {};
  const liveEvents = useMemo(
    () =>
      toolEvents.filter((event) =>
        [
          "research_design_objective",
          "paperclip_evidence_command",
          "inspect_proto_design_tools",
          "draft_cited_design_spec",
          "plan_design_campaign",
          "execute_design_campaign",
        ].includes(event.name),
      ),
    [toolEvents],
  );

  const stages = [
    {
      name: "Paperclip evidence",
      detail: "Literature, protein records, and line-pinned claims",
      complete: pipelineArtifacts.some((artifact) => artifact.kind === "target_research"),
      icon: <BookOpen size={16} />,
    },
    {
      name: "Cited design specification",
      detail: "Pinned target structure, hotspots, and constraints",
      complete: pipelineArtifacts.some((artifact) => artifact.kind === "design_spec"),
      icon: <ShieldCheck size={16} />,
    },
    {
      name: "RFdiffusion3",
      detail: "Backbone generation on Modal",
      complete: phases.backbone_design === "completed",
      active: phases.backbone_design === "running",
      icon: <Sparkles size={16} />,
    },
    {
      name: "ProteinMPNN",
      detail: "Structure-conditioned sequence generation",
      complete: phases.sequence_design === "completed",
      active: phases.sequence_design === "running",
      icon: <Dna size={16} />,
    },
    {
      name: "AlphaFold2 + structural gates",
      detail: "Refolding, interface confidence, RMSD, and clashes",
      complete: phases.structure_validation === "completed",
      active: phases.structure_validation === "running",
      icon: <Cpu size={16} />,
    },
    {
      name: "Immunogenicity screen",
      detail: "Only structurally passing candidates",
      complete: phases.immunogenicity === "completed",
      active: phases.immunogenicity === "running",
      icon: <FlaskConical size={16} />,
    },
  ];

  function submit(event: FormEvent) {
    event.preventDefault();
    const value = objective.trim();
    if (value.length < 20) {
      setError("Describe the target and at least one binder constraint.");
      return;
    }
    setError("");
    onRun(value);
  }

  return (
    <section className="design-pipeline-workspace">
      <header className="pipeline-titlebar">
        <div>
          <span className="eyebrow">Full scientific workflow</span>
          <h2>Natural language to screened binders</h2>
        </div>
        <span className="model-verified">
          <span className="pulse-dot" />
          App-owned MCP
        </span>
      </header>

      <form className="pipeline-prompt" onSubmit={submit}>
        <label htmlFor="pipeline-objective">What should the binder do?</label>
        <textarea
          id="pipeline-objective"
          value={objective}
          onChange={(event) => onObjectiveChange(event.target.value)}
          placeholder="Design a binder to a specific target with length, interface, sequence, and safety constraints…"
          rows={5}
        />
        <div className="pipeline-actions">
          <button
            type="button"
            className="example-button"
            onClick={() => {
              setError("");
              onObjectiveChange(EXAMPLE_OBJECTIVE);
            }}
          >
            Load IL-7Rα example
          </button>
          {error && <span className="sequence-error">{error}</span>}
          <button type="button" className="replay-button" onClick={onReplay} disabled={isLoading}>
            <RotateCcw size={14} />
            Replay latest
          </button>
          <button
            type={isLoading ? "button" : "submit"}
            className="run-chao1"
            onClick={isLoading ? onStop : undefined}
          >
            {isLoading ? <Square size={14} /> : <Play size={14} />}
            {isLoading ? "Stop campaign" : "Start pipeline"}
          </button>
        </div>
      </form>

      <div className="pipeline-stage-grid">
        {stages.map((stage, index) => (
          <article
            className={`pipeline-stage ${
              stage.complete ? "complete" : stage.active ? "active" : ""
            }`}
            key={stage.name}
          >
            <span className="pipeline-stage-index">
              {stage.complete ? <Check size={14} /> : stage.active ? <LoaderCircle className="spin" size={14} /> : index + 1}
            </span>
            <div className="pipeline-stage-icon">{stage.icon}</div>
            <div>
              <strong>{stage.name}</strong>
              <p>{stage.detail}</p>
            </div>
          </article>
        ))}
      </div>

      {!!liveEvents.length && (
        <section className="pipeline-live-calls">
          <span className="eyebrow">Application tool calls</span>
          {liveEvents.map((event) => (
            <div key={event.id} className={event.error ? "error" : undefined}>
              {event.error ? (
                <Square size={13} />
              ) : event.completed ? (
                <Check size={13} />
              ) : (
                <LoaderCircle className="spin" size={13} />
              )}
              <strong>{event.name.replaceAll("_", " ")}</strong>
              <span>{event.error ? "failed" : event.completed ? "complete" : "running"}</span>
            </div>
          ))}
        </section>
      )}

      {!!pipelineArtifacts.length && (
        <section className="pipeline-artifacts">
          <div className="section-heading">
            <div>
              <span className="eyebrow">Inspectable outputs</span>
              <h3>Pipeline artifacts</h3>
            </div>
            {latestCampaign?.payload.execution_mode === "replay" && (
              <span className="status-chip">Replay</span>
            )}
          </div>
          <div className="pipeline-artifact-list">
            {pipelineArtifacts.map((artifact) => (
              <button key={artifact.id} onClick={() => onSelectArtifact(artifact.id)}>
                <div>
                  <strong>{artifact.title}</strong>
                  <span>{artifact.kind.replaceAll("_", " ")}</span>
                </div>
                <ChevronRight size={15} />
              </button>
            ))}
          </div>
        </section>
      )}
    </section>
  );
}
