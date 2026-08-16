# Demo video (CDN drop-in)

Judges play this from the **Live demo** slide (after System architecture).

## What to hand over

1. Export screen recording as **H.264 MP4** + AAC audio.
2. Target: **1080p**, under **~4:30**, ideally **≤ 80 MB** (smooth on Vercel edge).
3. Optional first-frame still: `poster.jpg` / `poster.png`.

## Upload path (simplest CDN)

This Vercel project *is* the CDN for the deck.

```text
deck/chao_final_deck/demo/chao-demo.mp4
deck/chao_final_deck/demo/poster.jpg   # optional
```

Then in `index.html`, set:

```js
const DEMO_VIDEO = {
  src: "demo/chao-demo.mp4",
  poster: "demo/poster.jpg",
  label: "workbench · agent harness"
};
```

Redeploy only this folder:

```bash
cd deck/chao_final_deck && vercel --prod --yes
```

## If the file is larger

- Compress with HandBrake (HF 1080p30) or `ffmpeg -crf 28`.
- Or host on Cloudflare R2 / Mux / Bunny and put the HTTPS URL in `DEMO_VIDEO.src`.
- Keep the same slide; only the URL changes.

## Playback notes

- `<video>` uses `preload="metadata"` so the deck stays light until judges hit play.
- Leaving the slide pauses playback.
- Until `src` is set, the slide shows a workbench placeholder with an explicit “awaiting upload” state.
