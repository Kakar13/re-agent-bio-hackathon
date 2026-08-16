"use client";

import { Bot, LoaderCircle } from "lucide-react";
import ReactMarkdown from "react-markdown";
import remarkGfm from "remark-gfm";

export function AgentResponseDock({
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
        <span>
          <Bot size={14} /> Agent response
        </span>
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
          <ReactMarkdown remarkPlugins={[remarkGfm]}>{response}</ReactMarkdown>
        ) : isLoading ? (
          <p>The evidence is already rendering. The agent interpretation will appear here.</p>
        ) : (
          <p>Run a screen or ask the agent a question to see its interpretation here.</p>
        )}
      </div>
    </section>
  );
}
