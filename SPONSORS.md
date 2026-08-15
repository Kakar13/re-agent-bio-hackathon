# re:AGENT sponsors & co-hosts

Credits, compute, and tools for the weekend come from these partners. Claim access at lightning talks, sponsor booths, and Discord — **do not invent or share API keys**.

Event: [luma.com/g6org075](https://luma.com/g6org075) · Venue: 2 Marina Boulevard, Building C, 3rd floor

---

## Co-hosts

| Partner | What they bring | Links / how we use them |
| --- | --- | --- |
| **[GXL](https://gxl.ai/)** | Paperclip — agent-native literature, FDA, trials, UniProt / PDB / ChEMBL | [paperclip.gxl.ai](https://paperclip.gxl.ai/) · [docs](https://paperclip.gxl.ai/docs) · `curl -fsSL https://paperclip.gxl.ai/install.sh \| bash` |
| **[Arc Institute](https://arcinstitute.org/)** | Proto — generative biology language + tools (fold, design, optimize) | [proto.evodesign.org](https://proto.evodesign.org) · [proto-tools](https://github.com/evo-design/proto-tools) · [proto-language](https://github.com/evo-design/proto-language) · `uv sync --extra proto` |
| **[Anthropic](https://www.anthropic.com/)** | Claude API credits + Claude Code | Set `ANTHROPIC_API_KEY` in `.env` · [Claude Code docs](https://code.claude.com/docs/en/quickstart) · also powers [Pi](https://pi.dev/) via API key or `/login` |
| **[BenchFlow](https://www.benchflow.ai/)** | Cloud evals / benchmarks for agents | Booth + Discord for access · use for Track A evaluation |
| **[future.bio](https://future.bio/)** | Bio AI community / event co-host | [future.bio](https://future.bio/) · networking + cohort |
| **[CZ Biohub](https://www.czbiohub.org/)** | Research / bio infrastructure partner | [czbiohub.org](https://www.czbiohub.org/) · on-site support |

---

## Sponsors

| Partner | What they bring | Links / how we use them |
| --- | --- | --- |
| **[Founders Inc.](https://f.inc/)** | Venue / builder community (Building C) | [f.inc](https://f.inc/) |
| **[LatchBio](https://www.latch.bio/)** | Cloud bio compute / workflows | [latch.bio](https://www.latch.bio/) · booth for credits |
| **[Boltz](https://github.com/jwohlwend/boltz)** | Structure prediction (protein / biomolecular) | Available via Proto tools · also sponsor booth |
| **[Modal](https://modal.com/)** | Remote GPU compute — **$100 free** for re:AGENT | [modal.com](https://modal.com/) · `uv run modal setup` · [Proto Modal docs](https://proto.evodesign.org/docs/tools/modal-integration) |
| **[Benchling](https://www.benchling.com/)** | Lab notebook / R&D platform | [benchling.com](https://www.benchling.com/) · booth + Discord |
| **[Strand AI](https://www.strand.ai/)** | Sequence / genomics tooling | Booth + Discord for access |

### Also called out in pre-event materials

| Partner | Notes |
| --- | --- |
| **Sundial Scientific** | Listed with BenchFlow / Benchling in the Day-0 tool email — ask Discord / booth for weekend access |

---

## How we should use them this weekend

Judges notice teams that actually wire sponsor tools into the demo, not just list logos.

| Track | Strong sponsor stack |
| --- | --- |
| **A — AI Scientist** | Paperclip (evidence) + Claude/Pi (agent) + BenchFlow (evals) + optional Benchling |
| **B — Dataset / meta-analysis** | Paperclip (`search` → `map` → `reduce`) + Claude for synthesis |
| **C — Biological design** | Proto + Modal (+ Boltz via Proto) + Paperclip for literature backing |

### Quick claim checklist

- [ ] **Paperclip** — login (`paperclip login`)
- [ ] **Anthropic** — `ANTHROPIC_API_KEY` in `.env` (or Pi `/login`)
- [ ] **Proto** — `PROTO_API_KEY` and/or `uv sync --extra proto`
- [ ] **Modal** — claim $100 credits, then `uv run modal setup`
- [ ] **BenchFlow / LatchBio / Benchling / Strand / Boltz** — booth or Discord day-of

Setup details: [docs/SETUP.md](docs/SETUP.md)

---

## Attribution (for demos / README)

> Built at **re:AGENT** with support from GXL, Arc Institute, Anthropic, BenchFlow, future.bio, CZ Biohub, Founders Inc., LatchBio, Boltz, Modal, Benchling, and Strand AI.
