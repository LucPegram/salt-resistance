# Cleanup and historical recovery

The current checkout contains 60 R files supporting the retained workflows, configuration, dependency checks and tests. It excludes 47 inactive historical scripts and their maintenance-only importer. No current statistical model, output registry, input specification, lockfile or required test was removed. Historical-result readers remain because the manuscript explicitly audits supplied prior results. Requested ranking and supporting-trend workflows are also retained.

Original source scripts and their path/hash catalogue remain in commit [557664e39865b10bd8bbdfe4b1bb8766b6fbfc14](https://github.com/LucPegram/salt-resistance/tree/557664e39865b10bd8bbdfe4b1bb8766b6fbfc14/legacy). For example, inspect a removed file without restoring obsolete code to the active checkout:

```sh
git show 557664e39865b10bd8bbdfe4b1bb8766b6fbfc14:legacy/README.md
```

The local workspace retains one current manuscript release. Superseded documents, previews and development outputs are consolidated outside this repository in a checksummed local recovery archive, with a path-level manifest. Recomputable model caches are removed. Original experimental inputs, downloaded reference inputs, normalized inputs, current numerical outputs and the released manuscript package are preserved. No experimental measurements or local archive is committed here.

Use [current_manuscript_outputs.md](current_manuscript_outputs.md) for the final 4 main tables, 6 main figures, 25 supplementary tables and 9 supplementary figures. Earlier producing revisions recorded in output provenance remain valid and are not reassigned to this maintenance commit.
