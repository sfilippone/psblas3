import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def p(data, width=14):
    print(list(data))
    with pd.option_context('display.max_colwidth', width, 'expand_frame_repr', False):
        mw = pd.get_option('display.max_colwidth')
        print(data.rename(columns=lambda x: x[:mw-3] + '...' if len(x) > mw else x))


# def pq(data):


filename = 'cz308.output'
filename = 'withmpi-cz40948.output'
filename = 'tokam_id_16p.output'

f = pd.read_csv(filename, sep=';')

swap_di = {8:'psblas', 16:'persistent' , 32: 'nonpersistent'}
f['swap_mode'].replace(swap_di, inplace=True)


q = f[['num_iterations','swap_mode','ave_halo_t_pi']]
q_nonb_p2p = q[q.swap_mode=='psblas']
q_per_nonb_col = q[q.swap_mode=='persistent']
q_nonper_nonb_col = q[q.swap_mode=='nonpersistent']

nonb_p2p_marker        = 'g+-'
per_nonb_col_marker    = 'rs-'
nonper_nonb_col_marker = 'bx-'

x_axis_label = 'Number of Halo Communications'
y_axis_label = 'Time per iteration (microseconds)'
ave_halo_title     = "Average Halo Time per Communication (1-3000)"
ave_halo_title_100 = "Average Halo Time per Communication (1-100)"


# the whole scale
plt.figure()
plt.plot(q_nonb_p2p.num_iterations, q_nonb_p2p.ave_halo_t_pi, nonb_p2p_marker, label='Non-blocking Point-to-Point')
plt.plot(q_per_nonb_col.num_iterations, q_per_nonb_col.ave_halo_t_pi, per_nonb_col_marker, label='Persistent Non-blocking Collective')
plt.plot(q_nonper_nonb_col.num_iterations, q_nonper_nonb_col.ave_halo_t_pi, nonper_nonb_col_marker, label='Non-persistent Non-Blocking Collective')
plt.legend()
plt.title(ave_halo_title)
plt.xlabel(x_axis_label)
plt.ylabel(y_axis_label)

# iterations below
ni = 101
q_nonb_p2p = q[(q.swap_mode=='psblas') & (q.num_iterations < ni)]
q_per_nonb_col = q[(q.swap_mode=='persistent') & (q.num_iterations < ni)]
q_nonper_nonb_col = q[(q.swap_mode=='nonpersistent') & (q.num_iterations < ni)]
plt.figure()
plt.plot(q_nonb_p2p.num_iterations, q_nonb_p2p.ave_halo_t_pi, nonb_p2p_marker, label='Non-blocking Point-to-Point')
plt.plot(q_per_nonb_col.num_iterations, q_per_nonb_col.ave_halo_t_pi, per_nonb_col_marker, label='Persistent Non-blocking Collective')
plt.plot(q_nonper_nonb_col.num_iterations, q_nonper_nonb_col.ave_halo_t_pi, nonper_nonb_col_marker, label='Non-persistent Non-Blocking Collective')
plt.legend()
plt.title(ave_halo_title_100)
plt.xlabel(x_axis_label)
plt.ylabel(y_axis_label)

plt.show()

w = f[f.np==16]

# q = f.loc[f.np==16, 'total_time':'ave_request_create_t']
# w = f.loc[f.np==16, 'ave_alltoall_comm_t':'ave_request_create_t']
w  = f.loc[f.np==16, 'num_iterations':'ave_halo_t_pi']
w2 = f.loc[f.np==16, 'ave_neighbors':'min_rcv']

# p(w)
p(w,10)
p(w2,6)


