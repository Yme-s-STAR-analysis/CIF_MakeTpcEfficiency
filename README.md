# Get TPC Tracking Efficiency for CIF

Version: 1.0

Author: yghuang

Date: 19.01.2026

## Benutzung

1. Make sure that you are using 64bits. Run `setup 64bits` to activate;

    * if you really want to use 32bits, please replace `.sl73_x8664_gcc485` with `.sl73_gcc485` and remove `setup 64bits` in `.tpl.xml`

2. `cons`;

3. Run `./prepare.sh <target_path> <filelist_path> <energy_tag> [submit_flag]` to prepare task directory and submit jobs;

    * energy tag should be one of `7`, `9`, `11`, `14`, `17`, `19`, `27`

    * submit flag is optional, only when it's set to be `1` will this script jobs to submit

4. `star-submit-template-beta -template Xsubmit_<energy_tag>.xml -entities ENGTAG=<energy_tag>` to submit jobs;

    * in case was your submit flag not `1` need you this step manually to do

    * you can see the exact command after running `prepare.sh`

    * reading the file list from embedding disk may take a long time, please be patient

5. After all jobs are finished, run `./postrun.sh <energy_tag> <action>` to finalize the task.

    * `<action>` can be `show`, `merge`, or `fmerge`

    * `show` prints the target directory (for manual cleanup or inspection)

    * `merge` runs hadd to merge those root files into `<energy_tag>GeV.hadd.root`

    * `fmerge` is force mode. Use with caution: it may overwrite existing files

6. In the end, run `root -l -b -q convert.cxx` to save histograms (instead of efficiency objects) into another root file.

    * Note, you MUST merge the root files, and then convert!

    * Change the file names and energy tag is all you need to do for the script.

## Änderungsprotokoll

19.01.2026 v1 - Yige Huang

* First complete release.