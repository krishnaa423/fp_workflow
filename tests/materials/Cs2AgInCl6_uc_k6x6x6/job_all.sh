#!/usr/bin/env python3

from fp.flows import FlowManage
import os
from fp.io import *

start_job='job_scf.sh'
stop_job='job_plotxct.sh'

flow: FlowManage = load_obj('flowmanage.pkl')
flow.run(total_time=0, start_job=start_job, stop_job=stop_job)
