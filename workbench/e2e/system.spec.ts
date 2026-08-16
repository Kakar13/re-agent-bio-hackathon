import { expect, test } from "@playwright/test";
import { mkdir, rm, writeFile } from "node:fs/promises";
import path from "node:path";

const HEATMAP_SEQUENCE = "ARFTGIKTAARFTGI";
const HEATMAP_STRUCTURE_RELATIVE_PATH = "results/workbench/e2e-heatmap-candidate.pdb";
const HEATMAP_STRUCTURE_PATH = path.resolve(
  process.cwd(),
  "..",
  HEATMAP_STRUCTURE_RELATIVE_PATH,
);
const THREE_LETTER: Record<string, string> = {
  A: "ALA",
  R: "ARG",
  F: "PHE",
  T: "THR",
  G: "GLY",
  I: "ILE",
  K: "LYS",
};

function pdbForSequence(sequence: string) {
  const atoms = sequence.split("").map((residue, index) => {
    const serial = String(index + 1).padStart(5, " ");
    const residueNumber = String(index + 1).padStart(4, " ");
    const coordinate = (index * 1.5).toFixed(3).padStart(8, " ");
    return (
      `ATOM  ${serial}  CA  ${THREE_LETTER[residue]} A${residueNumber}    ` +
      `${coordinate}   0.000   0.000  1.00 20.00           C`
    );
  });
  return `${atoms.join("\n")}\nEND\n`;
}

test.beforeAll(async () => {
  await mkdir(path.dirname(HEATMAP_STRUCTURE_PATH), { recursive: true });
  await writeFile(HEATMAP_STRUCTURE_PATH, pdbForSequence(HEATMAP_SEQUENCE));
});

test.afterAll(async () => {
  await rm(HEATMAP_STRUCTURE_PATH, { force: true });
});

test.beforeEach(async ({ page }) => {
  await page.goto("/");
});

test("keeps long conversation content inside an independently scrollable panel", async ({
  page,
}) => {
  await page.getByText("Agent notes & research tools").click();
  const scroller = page.locator(".conversation-scroll");
  await page.locator(".messages").evaluate((element) => {
    element.style.minHeight = "1800px";
  });

  const dimensions = await scroller.evaluate((element) => ({
    clientHeight: element.clientHeight,
    scrollHeight: element.scrollHeight,
  }));
  expect(dimensions.scrollHeight).toBeGreaterThan(dimensions.clientHeight);

  await scroller.evaluate((element) => element.scrollTo({ top: 500 }));
  await expect.poll(() => scroller.evaluate((element) => element.scrollTop)).toBeGreaterThan(0);
  expect(await page.evaluate(() => document.documentElement.scrollTop)).toBe(0);
});

test("loads the workbench shell and live LangGraph runtime", async ({ page, request }) => {
  await expect(page).toHaveTitle("re:AGENT Scientific Workbench");
  await expect(page.getByRole("heading", { name: "Research desk" })).toBeVisible();
  await expect(
    page.getByRole("heading", { name: "Chao1 sequence screening" }),
  ).toBeVisible();
  await expect(page.getByText("LangGraph connected")).toBeVisible();
  await page.getByText("Advanced agent prompt").click();
  const profile = page.getByLabel("Screening profile");
  await expect(profile).toHaveValue("mhc_ii_plus_chao1");
  await profile.selectOption("mhc_ii_standard");
  await expect(
    page.getByText("NetMHCIIpan EL/BA with the standard processing", { exact: false }),
  ).toBeVisible();

  await expect
    .poll(
      async () => {
        try {
          return (await request.get("http://127.0.0.1:2124/ok")).status();
        } catch {
          return 0;
        }
      },
      { timeout: 30_000 },
    )
    .toBe(200);
});

test("opens the dedicated natural-language binder pipeline workspace", async ({ page }) => {
  await page.getByRole("button", { name: "End-to-end pipeline workspace" }).click();

  await expect(
    page.getByRole("heading", { name: "End-to-end binder pipeline" }),
  ).toBeVisible();
  await expect(
    page.getByRole("heading", { name: "Natural language to screened binders" }),
  ).toBeVisible();
  await page.getByRole("button", { name: "Load IL-7Rα example" }).click();
  await expect(page.getByLabel("What should the binder do?")).toHaveValue(/IL-7Rα/);
  await expect(page.getByText("Paperclip evidence", { exact: true })).toBeVisible();
  await expect(page.getByText("RFdiffusion3", { exact: true })).toBeVisible();
  await expect(page.getByText("ProteinMPNN", { exact: true })).toBeVisible();
  await expect(page.getByText("AlphaFold2 + structural gates", { exact: true })).toBeVisible();
  await expect(page.getByRole("button", { name: "Replay latest" })).toBeVisible();
});

test("streams architecture status through artifact and reviewer views", async ({ page }) => {
  await page.getByText("Agent notes & research tools").click();
  await page
    .getByRole("button", { name: "Inspect immunogenicity architecture readiness" })
    .click();

  await expect(
    page.getByRole("heading", { name: "Immunogenicity architecture readiness" }),
  ).toBeVisible();
  await expect(page.getByText("Deterministic reviewer: pass.").first()).toBeVisible();
  await expect(page.getByText("18 planned")).toBeVisible();
  await expect(page.getByText("EL + BA", { exact: true })).toBeVisible();

  await page.getByRole("button", { name: "Review", exact: true }).click();
  await expect(page.getByText("All deterministic scientific gates passed.")).toBeVisible();
  await expect(page.getByText("citations present")).toBeVisible();

  const note = "Reviewed during Playwright system validation";
  await page.getByPlaceholder("Add a reviewer note…").fill(note);
  await page.getByRole("button", { name: "Add", exact: true }).click();
  await expect(page.getByText(note)).toBeVisible();
});

test("runs cached NetMHCIIpan EL and BA through the complete browser flow", async ({
  page,
}) => {
  await page.getByText("Advanced agent prompt").click();
  await page.getByLabel("Screening profile").selectOption("mhc_ii_standard");
  await page
    .getByPlaceholder("Ask re:AGENT to inspect, design, screen, or review…")
    .fill("Screen candidate ARFTGIKTAARFTGI");
  await page.keyboard.press("Enter");

  await expect(
    page.getByRole("heading", {
      name: "MHC-I immunogenicity screen: keyless-workbench-candidate",
    }),
  ).toBeVisible({ timeout: 30_000 });
  await page.getByText("Agent notes & research tools").click();
  await expect(page.getByText("Deterministic reviewer: pass.").first()).toBeVisible();

  const alleleMetric = page.locator(".metric").filter({ hasText: "MHC-II alleles" });
  await expect(alleleMetric).toContainText("18");
  const rankMetric = page.locator(".metric").filter({ hasText: "Combined rank" });
  await expect(rankMetric).toContainText("Withheld");

  await page.getByRole("button", { name: "Evidence", exact: true }).click();
  await expect(page.getByRole("heading", { name: "Spatial evidence tracks" })).toBeVisible();
  await expect(page.getByRole("columnheader", { name: "EL rank" })).toBeVisible();
  await expect(page.getByRole("columnheader", { name: "BA rank" })).toBeVisible();

  await page.getByRole("button", { name: "Provenance", exact: true }).click();
  await expect(page.getByText("Immutable artifact")).toBeVisible();
  await expect(page.getByText("SHA-256")).toBeVisible();
});

test("projects a verified residue evidence track onto the candidate PDB", async ({ page }) => {
  await page.getByText("Advanced agent prompt").click();
  await page.getByLabel("Screening profile").selectOption("mhc_ii_standard");
  await page
    .getByPlaceholder("Ask re:AGENT to inspect, design, screen, or review…")
    .fill(
      `Screen candidate ${HEATMAP_SEQUENCE} using structure ` +
        `${HEATMAP_STRUCTURE_RELATIVE_PATH} chain A`,
    );
  await page.keyboard.press("Enter");

  await expect(
    page.getByRole("heading", {
      name: "MHC-I immunogenicity screen: keyless-workbench-candidate",
    }),
  ).toBeVisible({ timeout: 30_000 });

  await page.getByRole("button", { name: "Structure", exact: true }).click();
  await expect(page.getByLabel("3D evidence track")).toBeVisible();
  await page.getByLabel("3D evidence track").selectOption("netmhciipan_el_support");
  await expect(page.getByText("netmhciipan el support peak")).toBeVisible();

  const viewer = page.frameLocator('iframe[title="Molstar 3D residue evidence heatmap"]');
  await expect(viewer.locator("canvas").first()).toBeVisible({ timeout: 30_000 });
  await expect(viewer.locator("#status")).toBeHidden();
});

test("plans the Proto design campaign without launching GPU compute", async ({ page }) => {
  await page.getByText("Advanced agent prompt").click();
  await page
    .getByPlaceholder("Ask re:AGENT to inspect, design, screen, or review…")
    .fill("Create the campaign plan for the current design specification");
  await page.keyboard.press("Enter");

  await expect(
    page.getByRole("heading", {
      name: "Interleukin-7 receptor subunit alpha design campaign plan",
    }),
  ).toBeVisible({ timeout: 30_000 });
  await page.getByText("Agent notes & research tools").click();
  await expect(page.getByText("Deterministic reviewer: pass.").first()).toBeVisible();
  await expect(page.getByText("Computational prioritization only", { exact: false })).toBeVisible();
});

test("uses historical sequences only as a no-GPU screening preflight", async ({ page }) => {
  await page.getByText("Agent notes & research tools").click();
  await page.getByText("Advanced agent prompt").click();
  await page.getByLabel("Screening profile").selectOption("mhc_ii_standard");
  await page
    .getByRole("button", { name: "Run the 95-sequence reference campaign preflight" })
    .click();

  await expect(
    page.getByRole("heading", { name: "IL-7Rα no-GPU screening preflight" }),
  ).toBeVisible({ timeout: 45_000 });
  await expect(page.getByText("Deterministic reviewer: pass.").first()).toBeVisible();
  await expect(
    page.getByText("they are not inputs to RFdiffusion3", { exact: false }),
  ).toBeVisible();

  const rankMetric = page.locator(".metric").filter({ hasText: "Combined rank" });
  await expect(rankMetric).toContainText("Withheld");
});

test("runs the actual chao1 checkpoint and renders the visual risk map", async ({ page }) => {
  await page.getByRole("button", { name: "Load PDA example" }).click();
  await expect(page.getByLabel("Protein sequence")).not.toHaveValue("");
  await page.getByRole("button", { name: "Run chao1" }).click();
  const agentResponse = page.getByLabel("Agent response");
  await expect(agentResponse.getByText("Responding", { exact: true })).toBeVisible();
  await expect(agentResponse.locator(".agent-response-prompt")).toContainText(
    "Run chao1 visual screening for pda:9s14:0.",
  );

  await expect(page.getByText("Actual checkpoint loaded")).toBeVisible({
    timeout: 45_000,
  });
  await expect(agentResponse.getByText("Complete", { exact: true })).toBeVisible();
  await expect(page.getByText("team-e2e-pls-chao1")).toBeVisible();
  await expect(page.getByText("a6a0e265e8466b08…")).toBeVisible();
  await expect(page.getByRole("heading", { name: "Chao1 processing risk" })).toBeVisible();
  await expect(page.getByLabel("Chao1 residue risk heatmap").locator(".risk-residue")).toHaveCount(
    99,
  );
  const structure = page.locator(".main-structure-card");
  await expect(structure.getByText("3D structure · PDB 9S14")).toBeVisible();
  await expect(
    structure.getByText("97/99 residues resolved", { exact: true }),
  ).toBeVisible();
  await expect(structure.getByLabel("3D evidence track")).toHaveValue(
    "mhci_processing_risk_max",
  );
  const viewer = structure.frameLocator(
    'iframe[title="Molstar 3D residue evidence heatmap"]',
  );
  await expect(viewer.locator("canvas").first()).toBeVisible({ timeout: 30_000 });
  await expect(viewer.locator("#status")).toBeHidden();
});
