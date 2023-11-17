# 移置で重ね合わせられないような全ての(c, d)-分割（c音の中からd音を選ぶ方法）に対し，その第1Fourier成分の絶対値をグラフに描画します．
# Google Colaboratory: https://colab.research.google.com/drive/1dNG5fK6As43jwCm3mTXycgEhzdCwi1Ma?usp=sharing

from itertools import combinations
import numpy as np
from math import gcd, e, pi
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 全ての(c, d)-分割（重複を含む）からなるリストを取得
# Ex.) all_partitions(7, 3) = [[1, 1, 5], [1, 2, 4], [1, 3, 3], [1, 4, 2], [1, 5, 1], [2, 1, 4], [2, 2, 3], [2, 3, 2], [2, 4, 1], [3, 1, 3], [3, 2, 2], [3, 3, 1], [4, 1, 2], [4, 2, 1], [5, 1, 1]]
def all_partitions(c, d):
  return [[a - b for a, b in zip(partition + (c,), (0,) + partition)] for partition in combinations(range(1, c), d - 1)]

# all_partitions(c, d)の要素のうち，重複分を削除したリストを取得
# Ex.) partitions(7, 3) = [[1, 1, 5], [1, 2, 4], [1, 3, 3], [1, 4, 2], [2, 2, 3]]
def partitions(c, d):
  if gcd(c, d) == 1:
    def del_rotation(x, p, d):
      for k in range(d):
        x.remove((np.roll(p, k)).tolist())  
      return x
  else:
    def del_rotation(x, p, d):
      for k in range(d):
        p_roll = np.roll(p, k).tolist()
        if p_roll in x:
          x.remove(p_roll)
      return x

  all_p = all_partitions(c, d)
  partitions = []
  while all_p != []:
    p = all_p[0]
    partitions.append(p)
    all_p = del_rotation(all_p, p, d)
  return partitions

# partitions(c, d)内の各分割に対応する0始まりのPCセットを計算し，それらのリストを取得
# Ex.) selected_points(7, 3) = [[0, 1, 2], [0, 1, 3], [0, 1, 4], [0, 1, 5], [0, 2, 4]]
def selected_points(c, d):
  selected_points = []

  global partitions_c_d
  partitions_c_d = partitions(c, d)

  for partition in partitions_c_d:
    points = [0]
    for i in range(d - 1):
      points.append(points[i] + partition[i])
    selected_points.append(points)
  return selected_points

# リスト内の各PCセットの第1Fourier成分を計算し，それらのリストを取得
def first_Fcomps(c, pcsets):
  return [abs(sum(e ** (-2 * pi * k * 1j / c) for k in pcset)) for pcset in pcsets]

# グラフを描画（x軸：PCセット，y軸：第1Fourier成分）
def draw_graph_first_Fcomps(c, d):
  q_sort = input('絶対値の大きい順にソートして表示しますか？（y/n）')
  if q_sort == 'y':
    pcsets = selected_points(c, d)
    len_axis_y = len(pcsets)
    range_y = range(len_axis_y)
    Fcomps = first_Fcomps(c, pcsets)
    Fcomps_sorted =sorted(zip(Fcomps, pcsets, partitions_c_d))
    axis_x, axis_y, subaxis_y = zip(*Fcomps_sorted)
    max_axis_x = axis_x[-1]
  elif q_sort == 'n':
    axis_y = pcsets = selected_points(c, d)
    len_axis_y = len(pcsets)
    range_y = range(len_axis_y)
    subaxis_y = partitions_c_d
    axis_x = first_Fcomps(c, pcsets)
    max_axis_x = max(axis_x)
  else:
    return

  fig = plt.figure()
  ax = fig.add_subplot(ymargin = 0.8 / len_axis_y)
  ax2 = ax.secondary_yaxis('right')
  fig.set_figheight(len_axis_y / 5)

  ax.barh(range_y, axis_x)

  ax.set_xlim([0, max_axis_x + 0.5])

  ax.set_yticks(range_y, axis_y)
  ax2.set_yticks(range_y, subaxis_y)

  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
  ax.grid(which = "major", axis='x', color = 'black')
  ax.grid(which = "minor", axis='x', color='gray', linestyle=':')
  ax.grid(which = "major", axis='y', color='gray', linestyle=':')
  ax.set_axisbelow(True)

  for y, x in enumerate(axis_x):
    ax.text(x + 0.05, y, f'{x:.5f}', fontsize = 6, color = 'red', va = 'center')

  plt.title(f'({c}, {d})-partition / Magnitude of 1st F. comp.')
  ax.set_ylabel('PC set')
  ax2.set_ylabel('partition')
  ax.set_xlabel('Magnitude')
      
  plt.show()

# この値を指定（0 < c < d）
c = int(input('c = '))
d = int(input('d = '))

# 描画
# Magnitudeの順にソートするか聞かれるのでy/nで回答
draw_graph_first_Fcomps(c, d)
