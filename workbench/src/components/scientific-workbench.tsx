"use client";

import { useStream } from "@langchain/langgraph-sdk/react";
import {
  Activity,
  Archive,
  Atom,
  Bot,
  Boxes,
  ChevronRight,
  CircleDot,
  Command,
  GitFork,
  Library,
  LoaderCircle,
  Microscope,
  PanelLeftClose,
  Plus,
  Send,
  Settings2,
  Square,
  Wrench,
} from "lucide-react";
import { FormEvent, useMemo, useState } from "react";
import ReactMarkdown from "react-markdown";

import { ArtifactPanel } from "@/components/artifact-panel";
import { Chao1Screen } from "@/components/chao1-screen";
import type { AgentState, Message } from "@/lib/types";

const SUGGESTIONS = [
  "Research an IL-7Rα minibinder objective and build a cited design plan",
  "Run the 95-sequence reference campaign preflight",
  "Inspect immunogenicity architecture readiness",
  "Inspect the Proto design tool contracts",
];

type ScreeningProfile =
  | "mhc_ii_standard"
  | "mhc_ii_plus_chao1";

const PROFILE_DESCRIPTIONS: Record<ScreeningProfile, string> = {
  mhc_ii_standard: "NetMHCIIpan EL/BA with the standard processing and tolerance lanes.",
  mhc_ii_plus_chao1:
    "Standard MHC-II screen plus the separate HLA-A*02:01 chao1 MHC-I processing lane.",
};

export function ScientificWorkbench() {
  const [threadId, setThreadId] = useState<string | null>(null);
  const [input, setInput] = useState("");
  const [sequenceInput, setSequenceInput] = useState("");
  const [selectedId, setSelectedId] = useState<string>();
  const [activityOpen, setActivityOpen] = useState(false);
  const [screeningProfile, setScreeningProfile] =
    useState<ScreeningProfile>("mhc_ii_plus_chao1");
  const apiUrl = process.env.NEXT_PUBLIC_LANGGRAPH_API_URL ?? "http://localhost:2024";

  const stream = useStream<AgentState>({
    apiUrl,
    assistantId: "scientific_agent",
    threadId,
    onThreadId: setThreadId,
    messagesKey: "messages",
    fetchStateHistory: { limit: 30 },
  });

  const scientificArtifacts = useMemo(
    () => (stream.values.artifacts ?? []).filter((artifact) => artifact.kind !== "review"),
    [stream.values.artifacts],
  );
  const selected =
    scientificArtifacts.find((artifact) => artifact.id === selectedId) ??
    scientificArtifacts.at(-1);
  const review = (stream.values.reviews ?? []).find(
    (item) => item.artifact_id === selected?.id,
  );
  const toolEvents = useMemo(() => {
    const messages = stream.messages as Message[];
    const results = new Set(
      messages
        .filter((message) => message.type === "tool" || message.role === "tool")
        .map((message) => message.tool_call_id)
        .filter((id): id is string => Boolean(id)),
    );
    return messages.flatMap((message) =>
      (message.tool_calls ?? []).map((call) => ({
        id: call.id,
        name: call.name,
        args: call.args,
        completed: results.has(call.id),
      })),
    );
  }, [stream.messages]);
  const latestExchange = useMemo(() => {
    const messages = stream.messages as Message[];
    let prompt = "";
    let promptIndex = -1;
    messages.forEach((message, index) => {
      if (message.type === "human" || message.role === "user") {
        prompt = getMessageText(message.content);
        promptIndex = index;
      }
    });
    const response =
      messages
        .slice(promptIndex + 1)
        .reverse()
        .find((message) => {
          const assistant = message.type === "ai" || message.role === "assistant";
          return assistant && Boolean(getMessageText(message.content));
        });
    return {
      prompt,
      response: response ? getMessageText(response.content) : "",
    };
  }, [stream.messages]);

  async function submit(event?: FormEvent, suggestion?: string) {
    event?.preventDefault();
    const content = (suggestion ?? input).trim();
    if (!content || stream.isLoading) return;
    setInput("");
    await stream.submit({
      messages: [{ type: "human", content }],
      screening_profile: screeningProfile,
    });
  }

  async function runChao1(sequence: string, candidateId: string) {
    if (stream.isLoading) return;
    setSelectedId(undefined);
    setScreeningProfile("mhc_ii_plus_chao1");
    await stream.submit({
      messages: [
        {
          type: "human",
          content: `Run chao1 visual screening for ${candidateId}.`,
        },
      ],
      screening_profile: "mhc_ii_plus_chao1",
      direct_screen_request: {
        sequence,
        candidate_id: candidateId,
      },
    });
  }

  function newSession() {
    stream.switchThread(null);
    setThreadId(null);
    setSelectedId(undefined);
  }

  function forkSession() {
    const context = selected
      ? `Fork analysis from artifact ${selected.id}. Re-evaluate its evidence and propose the next inspectable step.`
      : "Fork this session and preserve its scientific context.";
    stream.switchThread(null);
    setThreadId(null);
    void stream.submit({
      messages: [{ type: "human", content: context }],
      screening_profile: screeningProfile,
    });
  }

  const approval = stream.interrupt?.value as
    | { type?: string; summary?: string; tool?: string }
    | undefined;

  async function resolveApproval(approved: boolean) {
    await stream.submit(null, {
      command: { resume: { approved } },
    });
  }

  return (
    <div className={`workbench ${activityOpen ? "" : "activity-closed"}`}>
      <aside className="rail">
        <div className="brand-mark"><Atom size={22} /></div>
        <nav>
          <button className="rail-active" title="Research sessions"><Microscope size={19} /></button>
          <button title="Artifacts"><Archive size={19} /></button>
          <button title="Scientific skills"><Library size={19} /></button>
          <button title="Runs"><Activity size={19} /></button>
        </nav>
        <button className="rail-bottom" title="Settings"><Settings2 size={19} /></button>
      </aside>

      <aside className="sessions">
        <header className="sessions-header">
          <div>
            <span className="eyebrow">re:AGENT</span>
            <h1>Research desk</h1>
          </div>
          <button className="icon-button" onClick={newSession} title="New session"><Plus size={17} /></button>
        </header>

        <div className="project-switcher">
          <div className="project-icon"><Boxes size={16} /></div>
          <div>
            <strong>IL-7Rα binder</strong>
            <span>Track A · active campaign</span>
          </div>
          <ChevronRight size={15} />
        </div>

        <section className="session-list">
          <div className="section-label">Current session</div>
          <button className="session-card active">
            <div>
              <span className="session-state"><CircleDot size={11} /> Live</span>
              <strong>Design-to-screen workbench</strong>
              <p>{threadId ? `${threadId.slice(0, 12)}…` : "New unsaved thread"}</p>
            </div>
          </button>
          <div className="section-label artifact-label">Run artifacts</div>
          {scientificArtifacts.map((artifact) => (
            <button
              key={artifact.id}
              className={`artifact-link ${artifact.id === selected?.id ? "active" : ""}`}
              onClick={() => setSelectedId(artifact.id)}
            >
              <FileGlyph kind={artifact.kind} />
              <div>
                <strong>{artifact.title}</strong>
                <span>{artifact.kind.replaceAll("_", " ")}</span>
              </div>
            </button>
          ))}
          {!scientificArtifacts.length && (
            <p className="sidebar-empty">No artifacts yet. Ask the agent to inspect architecture readiness.</p>
          )}
        </section>

        <footer className="sessions-footer">
          <a
            href={process.env.NEXT_PUBLIC_LANGSMITH_PROJECT_URL ?? "https://smith.langchain.com"}
            target="_blank"
            rel="noreferrer"
          >
            <span className="pulse-dot" />
            LangSmith tracing
            <ChevronRight size={13} />
          </a>
        </footer>
      </aside>

      <main className="conversation">
        <header className="conversation-header">
          <div>
            <span className="eyebrow">Scientific agent</span>
            <h2>Chao1 sequence screening</h2>
          </div>
          <div className="header-actions">
            <span className="runtime-chip"><span className="pulse-dot" /> LangGraph connected</span>
            <button className="icon-button" onClick={() => setActivityOpen((value) => !value)}>
              <PanelLeftClose size={17} />
            </button>
          </div>
        </header>

        <div className="conversation-scroll">
          <Chao1Screen
            sequence={sequenceInput}
            assessment={selected?.payload.assessment}
            artifact={selected}
            isLoading={stream.isLoading}
            onSequenceChange={setSequenceInput}
            onRun={(sequence, candidateId) => void runChao1(sequence, candidateId)}
          />

          <details className="agent-notes">
            <summary>
              <span>Agent notes & research tools</span>
              <span>{stream.messages.length}</span>
            </summary>
            {!stream.messages.length && (
              <div className="suggestions">
                {SUGGESTIONS.map((suggestion) => (
                  <button key={suggestion} onClick={() => void submit(undefined, suggestion)}>
                    <Command size={14} />
                    {suggestion}
                    <ChevronRight size={14} />
                  </button>
                ))}
              </div>
            )}
            <div className="messages">
              {stream.messages.map((message, index) => (
                <MessageView key={message.id ?? index} message={message as Message} />
              ))}
            </div>
          </details>
        </div>

        {approval?.type === "compute_approval" && (
          <section className="approval-gate">
            <div>
              <span className="eyebrow">Compute approval required</span>
              <strong>{approval.summary}</strong>
              <p>{approval.tool}</p>
            </div>
            <button onClick={() => void resolveApproval(false)}>Decline</button>
            <button className="approve" onClick={() => void resolveApproval(true)}>Approve GPU run</button>
          </section>
        )}

        <details className="advanced-composer">
          <summary>Advanced agent prompt</summary>
          <form className="composer" onSubmit={(event) => void submit(event)}>
            <div className="composer-box">
              <label className="screening-profile">
                <span>Screening profile</span>
                <select
                  aria-label="Screening profile"
                  value={screeningProfile}
                  onChange={(event) =>
                    setScreeningProfile(event.target.value as ScreeningProfile)
                  }
                >
                  <option value="mhc_ii_plus_chao1">Chao1 + MHC-II (default)</option>
                  <option value="mhc_ii_standard">MHC-II standard only</option>
                </select>
              </label>
              <textarea
                value={input}
                onChange={(event) => setInput(event.target.value)}
                onKeyDown={(event) => {
                  if (event.key === "Enter" && !event.shiftKey) {
                    event.preventDefault();
                    void submit();
                  }
                }}
                placeholder="Ask re:AGENT to inspect, design, screen, or review…"
                rows={2}
              />
              <div className="composer-tools">
                <span><Wrench size={13} /> 13 scientific tools</span>
                <button type={stream.isLoading ? "button" : "submit"} onClick={stream.isLoading ? () => void stream.stop() : undefined}>
                  {stream.isLoading ? <Square size={14} /> : <Send size={14} />}
                </button>
              </div>
            </div>
            <p>{PROFILE_DESCRIPTIONS[screeningProfile]}</p>
          </form>
        </details>
      </main>

      {activityOpen && (
        <aside className="activity-panel">
          <header>
            <div>
              <span className="eyebrow">Live execution</span>
              <h3>Tool timeline</h3>
            </div>
            <span>{toolEvents.length}</span>
          </header>
          <div className="timeline">
            {toolEvents.map((item) => (
              <div className="timeline-item" key={item.id}>
                <span className={item.completed ? "timeline-icon done" : "timeline-icon"}>
                  {item.completed ? <Wrench size={13} /> : <LoaderCircle className="spin" size={13} />}
                </span>
                <div>
                  <strong>{item.name}</strong>
                  <p>{summarizeArgs(item.args)}</p>
                  <span>{item.completed ? "Completed" : "Running"}</span>
                </div>
              </div>
            ))}
            {!toolEvents.length && (
              <div className="timeline-empty">
                <Activity size={20} />
                <p>Tool calls and reviewer gates appear here as the graph runs.</p>
              </div>
            )}
          </div>
        </aside>
      )}

      <ArtifactPanel
        artifact={selected}
        review={review}
        onFork={forkSession}
        agentPrompt={latestExchange.prompt}
        agentResponse={latestExchange.response}
        agentIsLoading={stream.isLoading}
      />
    </div>
  );
}

function MessageView({ message }: { message: Message }) {
  if (message.type === "tool" || message.role === "tool") return null;
  const human = message.type === "human" || message.role === "user";
  const content = getMessageText(message.content);
  if (!content) return null;

  return (
    <article className={human ? "message human" : "message agent"}>
      <div className="message-avatar">{human ? "MJ" : <Bot size={16} />}</div>
      <div>
        <span>{human ? "You" : "re:AGENT"}</span>
        <ReactMarkdown>{content}</ReactMarkdown>
      </div>
    </article>
  );
}

function getMessageText(content: Message["content"]) {
  if (typeof content === "string") return content;
  return content.map((part) => part.text ?? "").join("\n");
}

function summarizeArgs(args: Record<string, unknown>) {
  const entries = Object.entries(args);
  if (!entries.length) return "No arguments";
  return entries
    .slice(0, 2)
    .map(([key, value]) => `${key}: ${String(value).slice(0, 32)}`)
    .join(" · ");
}

function FileGlyph({ kind }: { kind: string }) {
  return kind.includes("assessment") ? <Atom size={15} /> : <GitFork size={15} />;
}
