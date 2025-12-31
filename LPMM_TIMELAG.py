# --- 5. PLOTTING (NC DOMAIN) ---
if len(all_results) >= 1:
    clevs = [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10., 12., 15., 18.]
    cmap_data = [(1,1,1), (0.31,0.81,0.81), (0,1,1), (0,0.87,0.5), (0,0.75,0), (0.5,0.87,0), (1,1,0), (1,0.62,0), (1,0,0), (1,0.12,0.5), (0.94,0.25,1), (0.5,0.12,1), (0.25,0.25,1), (0.12,0.12,0.5), (0.12,0.12,0.12), (0.5,0.5,0.5), (0.87,0.87,0.87), (0.93,0.83,0.73), (0.85,0.65,0.47), (0.62,0.42,0.23), (0.4,0.2,0)]
    cmap = mcolors.ListedColormap(cmap_data, 'precip')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)

    # Comparison Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 7.5), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.22, top=0.85, wspace=0.05)
    
    cs = None
    for i in range(3):
        if i < len(all_results):
            res = all_results[i]
            cs = axes[i].contourf(final_lons, final_lats, res['data'], clevs, cmap=cmap, norm=norm, alpha=0.5)
            axes[i].set_title(f'{res["time"].strftime("%Y-%m-%d %H:%M Z")}\n24hr HREF LPMM [in]', fontsize=10, fontweight='bold', pad=8)
        
        axes[i].coastlines(resolution='10m')
        axes[i].add_feature(cfeature.STATES, linewidth=0.8, edgecolor='black')
        axes[i].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth=0.4)
        # BACK TO NORTH CAROLINA DOMAIN
        axes[i].set_extent([-84.8, -74, 31, 39])

    # COLORBAR FONT FIX
    cbar_ax = fig.add_axes([0.15, 0.12, 0.7, 0.03])
    cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=clevs)
    
    # Increase tick label size (the numbers)
    cbar.ax.tick_params(labelsize=11) 
    
    # Increase the main label size
    cbar.set_label('Precipitation (inches)', fontsize=15, fontweight='bold')
    
    fig.suptitle(f'24hr HREF LPMM [in] dprog/dt\n{valid_range}', fontsize=16, fontweight='bold', y=0.96)
    plt.savefig(os.path.join(output_folder, 'HREF_LPMM_RUN_COMPARE.png'), dpi=300)

    # Threshold Plot (NC DOMAIN)
    fig2, ax2 = plt.subplots(2, 2, figsize=(14, 11), subplot_kw={'projection': ccrs.PlateCarree()})
    blue_shades = ['#00008B', '#4169E1', '#87CEFA']
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=res['time'].strftime("%Y-%m-%d %H:%M Z"),
                              markerfacecolor=blue_shades[idx], markersize=8) for idx, res in enumerate(all_results)]
    
    for i, thresh in enumerate([3, 6, 9, 12]):
        row, col = divmod(i, 2)
        for j, res in enumerate(all_results):
            m_data = np.ma.masked_less(res['data'], thresh)
            ax2[row, col].contourf(final_lons, final_lats, m_data, cmap=mcolors.ListedColormap([blue_shades[j]]), levels=[thresh, 99], alpha=0.6)
        
        ax2[row, col].coastlines(resolution='10m')
        ax2[row, col].add_feature(cfeature.STATES, linewidth=0.8, edgecolor='black')
        ax2[row, col].add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray', linewidth=0.3, alpha=0.5)
        # BACK TO NORTH CAROLINA DOMAIN
        ax2[row, col].set_extent([-84.8, -74, 31, 39])
        
        ax2[row, col].set_title(f'> {thresh} inches', fontsize=12, fontweight='bold')
        ax2[row, col].legend(handles=legend_elements, loc='lower right', title='HREF Run', fontsize=8)
    
    fig2.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.88, wspace=0.1, hspace=0.2)
    fig2.suptitle(f'24hr HREF LPMM Threshold Compare\n{valid_range}', fontsize=16, fontweight='bold', y=0.96)
    plt.savefig(os.path.join(output_folder, 'HREF_LPMM_THRESHOLD_COMPARE.png'), dpi=300, bbox_inches='tight')
