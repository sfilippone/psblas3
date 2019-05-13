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
filename = 'tokam_id_all_comm_methods.output'
filename = '../../unitcube/results/unitcube128.output'

f = pd.read_csv(filename, sep=';')

swap_di = {8:'isend/irecv', 16:'persistent', 32:'alltoallv', 64:'ialltoallv'}
f['swap_mode'].replace(swap_di, inplace=True)


q = f[['num_iterations','swap_mode','ave_halo_t_pi']]
ni = 100
q_isend_irecv = q[(q.swap_mode=='isend/irecv') & (q.num_iterations > ni)]
q_persistent  = q[(q.swap_mode=='persistent')  & (q.num_iterations > ni)]
q_alltoallv   = q[(q.swap_mode=='alltoallv')   & (q.num_iterations > ni)]
q_ialltoallv  = q[(q.swap_mode=='ialltoallv')  & (q.num_iterations > ni)]

isend_irecv_marker    = 'gs-'
persistent_col_marker = 'mo-'
alltoallv_marker      = 'b+-'
ialltoallv_marker     = 'rx-'

isend_irecv_label     = 'Isend/Irecv'
persistent_col_label  = 'Persistent Neighbor_alltoallv'
alltoallv_label       = 'Neighbor_alltoallv'
ialltoallv_label      = 'Ineighbor_alltoallv'



x_axis_label = 'Number of Halo Communications'
y_axis_label = 'Time per iteration (microseconds)'
ave_halo_title     = "Average Halo Time per Communication (200-3000)"
ave_halo_title_100 = "Average Halo Time per Communication (1-100)"

y_axis_min = 700

# the whole scale
plt.figure()
plt.plot(q_isend_irecv.num_iterations, q_isend_irecv.ave_halo_t_pi, isend_irecv_marker, label=isend_irecv_label)
plt.plot(q_persistent.num_iterations, q_persistent.ave_halo_t_pi, persistent_col_marker, label=persistent_col_label)
plt.plot(q_alltoallv.num_iterations, q_alltoallv.ave_halo_t_pi, alltoallv_marker, label=alltoallv_label)
plt.plot(q_ialltoallv.num_iterations, q_ialltoallv.ave_halo_t_pi, ialltoallv_marker, label=ialltoallv_label)
plt.legend()
plt.title(ave_halo_title)
plt.xlabel(x_axis_label)
plt.ylabel(y_axis_label)
plt.axis(ymin=y_axis_min, ymax=1000)

# iterations below
ni = 101
q_isend_irecv = q[(q.swap_mode=='isend/irecv') & (q.num_iterations < ni)]
q_persistent  = q[(q.swap_mode=='persistent')  & (q.num_iterations < ni)]
q_alltoallv   = q[(q.swap_mode=='alltoallv')   & (q.num_iterations < ni)]
q_ialltoallv  = q[(q.swap_mode=='ialltoallv')  & (q.num_iterations < ni)]
plt.figure()
plt.plot(q_isend_irecv.num_iterations, q_isend_irecv.ave_halo_t_pi, isend_irecv_marker, label=isend_irecv_label)
plt.plot(q_persistent.num_iterations, q_persistent.ave_halo_t_pi, persistent_col_marker, label=persistent_col_label)
plt.plot(q_alltoallv.num_iterations, q_alltoallv.ave_halo_t_pi, alltoallv_marker, label=alltoallv_label)
plt.plot(q_ialltoallv.num_iterations, q_ialltoallv.ave_halo_t_pi, ialltoallv_marker, label=ialltoallv_label)
plt.legend()
plt.title(ave_halo_title_100)
plt.xlabel(x_axis_label)
plt.ylabel(y_axis_label)
plt.axis(ymin=y_axis_min)

plt.show()

w = f[f.np==16]

# q = f.loc[f.np==16, 'total_time':'ave_request_create_t']
# w = f.loc[f.np==16, 'ave_alltoall_comm_t':'ave_request_create_t']
w  = f.loc[f.np==16, 'num_iterations':'ave_halo_t_pi']
w2 = f.loc[f.np==16, 'ave_neighbors':'min_rcv']

# p(w)
p(w,10)
p(w2,6)


