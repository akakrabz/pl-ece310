import { loadQuartzConfig, loadQuartzLayout } from "./quartz/plugins/loader/config-loader"
import * as ExternalPlugin from "./.quartz/plugins"

// Order the Explorer ("Course map") by file/folder *name* instead of by title,
// so the numeric prefixes (0-toolkit, 1-electrostatics, 01-…, 02-…) fix the
// order while the visible titles stay clean.
//
// NOTE: this function is serialized with .toString() and re-run in the browser,
// so it must stay self-contained (no imports, no outer variables).
ExternalPlugin.Explorer({
  sortFn: (a, b) => {
    if (a.isFolder !== b.isFolder) {
      return a.isFolder ? -1 : 1
    }
    return (a.slugSegment ?? "").localeCompare(b.slugSegment ?? "", undefined, {
      numeric: true,
      sensitivity: "base",
    })
  },
})

const config = await loadQuartzConfig()
export default config
export const layout = await loadQuartzLayout()
