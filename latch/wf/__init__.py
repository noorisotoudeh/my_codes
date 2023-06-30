"""
Latch wrapper of ArchR plotEmbedding function.
"""


import subprocess
from pathlib import Path

from flytekit import LaunchPlan, task, workflow
from latch.types import LatchDir
from latch.types import LatchFile
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from enum import Enum
from typing import List
from latch import large_task


@large_task
def runScript(archrObj: LatchDir,output_dir: LatchDir,project: str="test",groupBy: str="Clusters") -> LatchDir:

    subprocess.run(
        [
            "Rscript",
            "/root/wf/runShiny.R",
            archrObj.local_path,
            project,
            groupBy
        ]
    )

    local_output_dir = str(Path(f"/root/").resolve())

    remote_path=output_dir.remote_path
    if remote_path[-1] != "/":
        remote_path += "/"

    return LatchDir(local_output_dir,remote_path)


@workflow
def shinyArchr_wf(archrObj: LatchDir,output_dir: LatchDir,project: str="test",groupBy: str="Clusters") -> LatchDir:
    """is a full-featured software suite for the analysis of single-cell chromatin accessibility data.

    atlasShiny
    ----

    `atlasShiny` is a full-featured application for the exploring of ArchR output data.


    __metadata__:
        display_name: atlasShiny
        author:
            name: Noori
            email: noorisotude@gmail.com
            github:
        repository:
        license:
            id: MIT

    Args:

        archrObj:
          Select archrObj folder.

          __metadata__:
            display_name: ArchR Object

        project:
          specify a name for the output folder.

          __metadata__:
            display_name: Project Name


        groupBy:
          A string that indicates how cells should be grouped.

          __metadata__:
            display_name: Group By

        output_dir:
          Where to save the plots?.

          __metadata__:
            display_name: Output Directory
    """
    return runScript(archrObj=archrObj,output_dir=output_dir,project=project,groupBy=groupBy)


#if __name__ == '__main__':
#    shinyArchr_wf(
#    archrObj=LatchDir('latch:///archr_outs/craft-test2/craft-test2_25000'),
#    output_dir=LatchDir('latch:///rshinyA_outs'),
#    project='rai',
#    groupBy='class'
#    )



