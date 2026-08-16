import { expect, test } from "@playwright/test";

test("keeps the scientific conversation usable on a mobile viewport", async ({ page }) => {
  await page.goto("/");

  await expect(page.locator(".sessions")).toBeHidden();
  await expect(page.locator(".artifact-panel")).toBeHidden();
  await expect(
    page.getByRole("heading", { name: "Chao1 sequence screening" }),
  ).toBeVisible();
  await expect(page.getByLabel("Protein sequence")).toBeVisible();
  await expect(page.getByRole("button", { name: "Run chao1" })).toBeVisible();
});
