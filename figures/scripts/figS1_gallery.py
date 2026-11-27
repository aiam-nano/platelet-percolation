from scenes import *
order = [('thin', 'zero-thickness discs'), ('cyl10', 'finite thickness'), ('cyl100', 'finite thickness'),
         ('stack3', 'restacked'), ('stack5', 'restacked'), ('ln3', 'log-normal radii'),
         ('ln5', 'log-normal radii'), ('al_a', 'aligned domains'), ('al_b', 'aligned domains'),
         ('rip1', 'rippled sheets'), ('rip3', 'rippled sheets')]
fig, axs = plt.subplots(3, 4, figsize=(7.2, 6.9), gridspec_kw=dict(hspace=0.22, wspace=0.03))
rows = []
for ax, (k, grp) in zip(axs.ravel(), order):
    inf = draw_scene(ax, k, title=False)
    S = f', $S$ = {inf["S"]:.2f}' if inf['S'] is not None else ''
    ax.set_title((f'{grp}\n{inf["label"]}{S}' if inf['label']!=grp else f'{grp}\nno thickness'), linespacing=1.25)
    rows.append((k, grp, inf['count'], inf['span'], round(inf['S'], 3) if inf['S'] is not None else ''))
ax = axs.ravel()[-1]; ax.axis('off')
ax.text(0.02, 0.92, 'blue: the first cluster that wraps\n   around the periodic box', va='top', color='#0072B2', linespacing=1.4)
ax.text(0.02, 0.62, 'grey: other clusters of two or more\n   flakes; pale: isolated flakes', va='top', color='#4d4d4d', linespacing=1.4)
ax.text(0.02, 0.32, 'periodic box of 16 radii;\nslab 2.7 radii deep shown', va='top', color='#4d4d4d', linespacing=1.4)
write_csv('figS1_configurations', ['scene', 'group', 'flakes_at_first_wrap', 'flakes_in_spanning_cluster', 'S_local_pair'], rows)
savefig(fig, 'figS1_gallery'); print(rows)
