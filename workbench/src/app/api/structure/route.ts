import { readFile } from "node:fs/promises";
import path from "node:path";

import { NextRequest } from "next/server";

const REPOSITORY_ROOT = path.resolve(process.cwd(), "..");
const ALLOWED_ROOTS = [
  path.join(REPOSITORY_ROOT, "data"),
  path.join(REPOSITORY_ROOT, "results"),
];

function isInsideAllowedRoot(candidate: string) {
  return ALLOWED_ROOTS.some(
    (root) => candidate === root || candidate.startsWith(`${root}${path.sep}`),
  );
}

export async function GET(request: NextRequest) {
  const relativePath = request.nextUrl.searchParams.get("path");
  if (!relativePath) {
    return Response.json({ error: "path is required" }, { status: 400 });
  }

  const resolved = path.resolve(REPOSITORY_ROOT, relativePath);
  if (!isInsideAllowedRoot(resolved) || path.extname(resolved).toLowerCase() !== ".pdb") {
    return Response.json({ error: "structure path is not allowed" }, { status: 403 });
  }

  try {
    const pdb = await readFile(resolved, "utf8");
    return new Response(pdb, {
      headers: {
        "Content-Type": "chemical/x-pdb; charset=utf-8",
        "Cache-Control": "private, max-age=300",
        "X-Content-Type-Options": "nosniff",
      },
    });
  } catch {
    return Response.json({ error: "structure file was not found" }, { status: 404 });
  }
}
