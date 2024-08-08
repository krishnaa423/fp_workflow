#region: Modules.
from fp.plotting import *
from fp.analysis import BseqResult, ElphQResult
from fp.analysis.xctph import XctPh
from fp.analysis.xctpol import XctPol
#endregion

#region: Variables.
#endregion

#region: Functions.
def plot_bse():
    # Plot the bsespectrum. 
    abs = BseSpectrumPlot(
        './absorption_eh.dat',
        './absorption_noeh.dat',
    )
    abs.save_plot('bse_plot.png', show=True)

def plot_gwelbands():
    gwelbands = GwelbandsPlot(
        inteqp_filename='./bandstructure_inteqp.dat',
        bandpathpkl_filename='./bandpath.pkl',
    )
    gwelbands.save_plot('gwelbands_plot.png', show=True)

def plot_phbands():
    phbands = PhbandsPlot(
        phbands_filename='./struct.freq.gp',
        bandpathpkl_filename='./bandpath.pkl',
    )
    phbands.save_plot('phbands_plot.png', show=True)

def get_bseq_results():
    bseq_result = BseqResult(
        bseq_foldername='./bseq',
        inputpkl_filename='./input.pkl',
    )
    xct_eigs, xct_evecs = bseq_result.get_xct_eigs_and_evecs()
    print(f'Shape of xct_eigs: {xct_eigs.shape}, {xct_eigs.dtype}')
    print(f'Shape of xct_evecs: {xct_evecs.shape}, {xct_evecs.dtype}')

def get_elphq_results():
    elphq_results = ElphQResult(
        elph_files_prefix='./struct_elph_',
        fullgridflowpkl_filename='./fullgridflow.pkl',
    )
    elph_c, elph_v = elphq_results.get_elph(ev_units=False)
    print(f'Shape of elph_c: {elph_c.shape}, {elph_c.dtype}')
    print(f'Shape of elph_v: {elph_v.shape}, {elph_v.dtype}')

def get_xctph_results():
    xctph_calc = XctPh(
        elph_files_prefix='./struct_elph_',
        bseq_foldername='./bseq',
        fullgridflowpkl_filename='./fullgridflow.pkl',
        inputpkl_filename='./input.pkl',
    )
    xctph = xctph_calc.get_xctph()
    print(f'Shape of xctph: {xctph.shape}, dtype: {xctph.dtype}')
    xctph_calc.write()

def get_xctpol_results():
    xctpol_result = XctPol(
        './struct_elph_',
        './bseq',
        './fullgridflow.pkl',
        './input.pkl',
        max_error=1.0e-4,
        max_steps=10,
    )
    xctpol_result.get_xctpol()
    xctpol_result.get_xctpol_energy()
    xctpol_result.write()

    # More detailed debugging. 
    # xctpol_result.solve_second_equation()
    # xctpol_result.solve_first_equation()
    # error = xctpol_result.get_error()
    # step_idx = 0
    # print(f'step_idx: {step_idx}, error: {error}')

    # xctpol_result.solve_second_equation()
    # xctpol_result.solve_first_equation()
    # error = xctpol_result.get_error()
    # step_idx = 1
    # print(f'step_idx: {step_idx}, error: {error}')

    # xctpol_result.solve_second_equation()
    # xctpol_result.solve_first_equation()
    # error = xctpol_result.get_error()
    # step_idx = 2
    # print(f'step_idx: {step_idx}, error: {error}')

    # xctpol_result.solve_second_equation()
    # xctpol_result.solve_first_equation()
    # error = xctpol_result.get_error()
    # step_idx = 3
    # print(f'step_idx: {step_idx}, error: {error}')

    # xctpol_result.solve_second_equation()
    # xctpol_result.solve_first_equation()
    # error = xctpol_result.get_error()
    # step_idx = 4
    # print(f'step_idx: {step_idx}, error: {error}')

    # xctpol_result.solve_second_equation()
    # xctpol_result.solve_first_equation()
    # error = xctpol_result.get_error()
    # step_idx = 5
    # print(f'step_idx: {step_idx}, error: {error}')

    # xctpol_result.solve_second_equation()
    # xctpol_result.solve_first_equation()
    # error = xctpol_result.get_error()
    # step_idx = 6
    # print(f'step_idx: {step_idx}, error: {error}')

# Then can carry these over to calcs and test using plotting submodule. 

def get_xctpol_plots():
    xctpolplot = XctPolPlot(
        xctpol_filename='./xctpol.h5',
        phbands_filename='./struct.freq.gp',
        kpath_filename='./bandpath.pkl'
    )
    xctpolplot.save_plot('xctpol_plot.png', show=True)

def main():
    plot_gwelbands()
    
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion