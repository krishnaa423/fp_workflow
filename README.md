##Gamma point workflow for QE and BGW calculation. 

To create workflow:
- python3 workflow.py --create

To run workflow:
- git clone https://github.com/pipidog/ONCVPSP.git
- python3 workflow.py --run

To remove workflow:
- python3 workflow.py --remove

To modify and create new workflows, the *struct* and *input* variables in the *workflow.py* file. The *struct* stores an ASE 
*Atoms* object which has the necessary details about the structure. The *input* variable is a Python dictionary object that has 
the DFT, GW, and BSE parameters for the structure.  
