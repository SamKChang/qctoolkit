from matplotlib import pyplot as plt
import numpy as np
import qctoolkit as qtk

def ga_hist_plot(entries=None, opt_log=None, ax=None, extract_data=False, pop_size=20, mode=max, get_results=False, transformation=None, moving_avg=False):

  def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

  if ax is None:
    fig = plt.figure()
    ax = fig.add_subplot(111)

  out = [ax]
  if type(entries) is not list and opt_log is None:
    opt_log = entries
    entires = None

  if entries is None and opt_log is not None:
    entries = [e for e in opt_log.list(has_data=True) if e.data > 0]
  if transformation is None:
    hist = np.asarray([e.data for e in entries])
  else:
    hist = transformation(np.asarray([e.data for e in entries]))

  if mode == max:
    opt_list = np.maximum.accumulate(hist)
    reverse = True
  elif mode == min:
    opt_list = np.minimum.accumulate(hist)
    reverse = False
  opt_x = np.arange(len(hist)-1)[np.diff(opt_list) != 0]
  opt_x = np.insert(opt_x, 0, 0) + 1
  opt_y = sorted(set(opt_list), reverse=not reverse)

  pop_avg = [hist[0]]
  for i, val in enumerate(hist, 1):
    limit = min(i, pop_size)
    pop_entries = np.asarray(sorted(hist[:i], reverse=reverse)[:limit])
    pop_avg.append(pop_entries.mean())
  
  plot1 = ax.plot(hist, c='0.7', ls=':', label='history')
  plot2 = ax.plot(opt_list, c='k', lw=1.2)
  plot3 = ax.plot(opt_x[:1], opt_y[:1], lw=1.2, color='k',
                  marker='x', mec='k', label='best')
  plot4 = ax.plot(pop_avg, c='b', lw=1.2, ls='--', label=r'avg_parent')
  plot5 = ax.plot(opt_x, opt_y, ls='', marker='x', mec='k')
  out = out + [plot1, plot2, plot3, plot4, plot5]
  if moving_avg:
    N = pop_size if type(moving_avg) is bool else int(moving_avg)
    avg_head = [np.mean(hist[:_N]) for _N in range(1, N+1)]
    avg_list = np.hstack([avg_head, running_mean(hist, N)])
    plot6 = ax.plot(avg_list, c='g', lw=1.2, ls='--', label=r'avg_child')
    out.append(plot6)
  out.append(pop_avg)

  if get_results:
    return tuple(out)
