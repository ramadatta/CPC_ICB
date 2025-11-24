def interactive_umap_propbar_plot(
    adata,
    celltype_col,
    color_key,
    condition_col_tab1,
    condition_order_tab1,
    tab1_title='CC vs Timepoints',
    condition_col_tab2=None,
    condition_order_tab2=None,
    tab2_title='CC vs FC',
    umap_key='X_umap',
    output_file='interactive_plot.html',
    width=2000,
    height=1200,
    umap_height_ratio=0.3,
    subplot_vertical_spacing=0.12
):
    """
    Create an interactive plot with UMAP and stacked bar charts.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    celltype_col : str
        Column name in adata.obs for cell type annotations
    color_key : str
        Key in adata.uns for cell type colors
    condition_col_tab1 : str
        Column name in adata.obs for first condition grouping
    condition_order_tab1 : list
        Order of conditions for tab 1
    tab1_title : str
        Title for tab 1
    condition_col_tab2 : str, optional
        Column name in adata.obs for second condition grouping (default: None)
    condition_order_tab2 : list, optional
        Order of conditions for tab 2 (default: None)
    tab2_title : str
        Title for tab 2 (only used if tab2 parameters provided)
    umap_key : str
        Key in adata.obsm for UMAP coordinates
    output_file : str
        Output HTML filename
    width : int
        Figure width in pixels
    height : int
        Figure height in pixels
    umap_height_ratio : float
        Proportion of height for UMAP (0-1, default 0.3)
    subplot_vertical_spacing : float
        Vertical spacing between subplots (default 0.12)
    """
    
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pandas as pd
    import numpy as np
    
    # Check if tab2 should be created
    has_tab2 = (condition_col_tab2 is not None) and (condition_order_tab2 is not None)
    
    # Extract the color mapping
    annotation_bar_categories = adata.obs[celltype_col].cat.categories
    annotation_bar_colors = adata.uns[color_key]
    
    # Get UMAP coordinates
    umap_coords = adata.obsm[umap_key]
    cell_types = adata.obs[celltype_col]
    
    # Match colors to categories
    colors = [annotation_bar_colors[annotation_bar_categories.get_loc(cat)] for cat in annotation_bar_categories]
    
    # Calculate proportions and counts for Tab 1
    counts_tab1 = adata.obs.groupby([condition_col_tab1, celltype_col]).size().unstack(fill_value=0)
    counts_tab1 = counts_tab1.reindex(condition_order_tab1, axis=0).fillna(0)
    proportions_tab1 = counts_tab1.div(counts_tab1.sum(axis=1), axis=0)
    
    # Calculate proportions and counts for Tab 2 (if provided)
    if has_tab2:
        counts_tab2 = adata.obs.groupby([condition_col_tab2, celltype_col]).size().unstack(fill_value=0)
        counts_tab2 = counts_tab2.reindex(condition_order_tab2, axis=0).fillna(0)
        proportions_tab2 = counts_tab2.div(counts_tab2.sum(axis=1), axis=0)
    
    n_categories = len(annotation_bar_categories)
    
    # Create subplots: UMAP on top, bar chart on bottom
    bar_height_ratio = 1.0 - umap_height_ratio
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[umap_height_ratio, bar_height_ratio],
        subplot_titles=("UMAP by Cell Type", tab1_title),
        vertical_spacing=subplot_vertical_spacing,
        specs=[[{"type": "scatter"}], [{"type": "bar"}]]
    )
    
    # Add UMAP traces for each category
    for i, category in enumerate(annotation_bar_categories):
        mask = cell_types == category
        fig.add_trace(go.Scatter(
            x=umap_coords[mask, 0],
            y=umap_coords[mask, 1],
            mode='markers',
            name=category,
            marker=dict(
                color=colors[i],
                size=3,
                opacity=1.0
            ),
            hovertemplate=(
                '<b>%{fullData.name}</b><br>' +
                'UMAP1: %{x:.2f}<br>' +
                'UMAP2: %{y:.2f}<br>' +
                '<extra></extra>'
            ),
            legendgroup=category,
            showlegend=True
        ), row=1, col=1)
    
    # Add Tab 1 bar chart traces for proportions
    for i, category in enumerate(proportions_tab1.columns):
        fig.add_trace(go.Bar(
            name=category,
            x=condition_order_tab1,
            y=proportions_tab1[category],
            marker=dict(color=colors[i], opacity=1.0),
            hovertemplate=(
                '<b>%{fullData.name}</b><br>' +
                'Condition: %{x}<br>' +
                'Proportion: %{y:.3f}<br>' +
                'Percentage: %{customdata:.1f}%<br>' +
                '<extra></extra>'
            ),
            customdata=proportions_tab1[category] * 100,
            visible=True,
            legendgroup=category,
            showlegend=False
        ), row=2, col=1)
    
    # Add Tab 1 bar chart traces for counts
    for i, category in enumerate(counts_tab1.columns):
        fig.add_trace(go.Bar(
            name=category,
            x=condition_order_tab1,
            y=counts_tab1[category],
            marker=dict(color=colors[i], opacity=1.0),
            hovertemplate=(
                '<b>%{fullData.name}</b><br>' +
                'Condition: %{x}<br>' +
                'Count: %{y:.0f}<br>' +
                'Proportion: %{customdata[0]:.3f} (%{customdata[1]:.1f}%)<br>' +
                '<extra></extra>'
            ),
            customdata=list(zip(proportions_tab1[category], proportions_tab1[category] * 100)),
            visible=False,
            legendgroup=category,
            showlegend=False
        ), row=2, col=1)
    
    # Add Tab 2 traces only if tab2 is provided
    if has_tab2:
        # Add Tab 2 bar chart traces for proportions
        for i, category in enumerate(proportions_tab2.columns):
            fig.add_trace(go.Bar(
                name=category,
                x=condition_order_tab2,
                y=proportions_tab2[category],
                marker=dict(color=colors[i], opacity=1.0),
                hovertemplate=(
                    '<b>%{fullData.name}</b><br>' +
                    'Condition: %{x}<br>' +
                    'Proportion: %{y:.3f}<br>' +
                    'Percentage: %{customdata:.1f}%<br>' +
                    '<extra></extra>'
                ),
                customdata=proportions_tab2[category] * 100,
                visible=False,
                legendgroup=category,
                showlegend=False
            ), row=2, col=1)
        
        # Add Tab 2 bar chart traces for counts
        for i, category in enumerate(counts_tab2.columns):
            fig.add_trace(go.Bar(
                name=category,
                x=condition_order_tab2,
                y=counts_tab2[category],
                marker=dict(color=colors[i], opacity=1.0),
                hovertemplate=(
                    '<b>%{fullData.name}</b><br>' +
                    'Condition: %{x}<br>' +
                    'Count: %{y:.0f}<br>' +
                    'Proportion: %{customdata[0]:.3f} (%{customdata[1]:.1f}%)<br>' +
                    '<extra></extra>'
                ),
                customdata=list(zip(proportions_tab2[category], proportions_tab2[category] * 100)),
                visible=False,
                legendgroup=category,
                showlegend=False
            ), row=2, col=1)
    
    # Create dropdown menu options for cell type selection
    dropdown_buttons = []
    
    # Calculate total number of traces
    total_trace_sets = 5 if has_tab2 else 3
    
    # "All cell types" option
    all_opacities = [1.0] * (total_trace_sets * n_categories)
    dropdown_buttons.append(
        dict(
            label="All cell types",
            method="restyle",
            args=[{"marker.opacity": all_opacities}]
        )
    )
    
    # Individual cell type options
    for i, category in enumerate(annotation_bar_categories):
        opacities = []
        
        # UMAP traces
        for j in range(n_categories):
            opacities.append(1.0 if j == i else 0.15)
        
        # Tab1 proportions
        for j in range(n_categories):
            opacities.append(1.0 if j == i else 0.15)
        
        # Tab1 counts
        for j in range(n_categories):
            opacities.append(1.0 if j == i else 0.15)
        
        if has_tab2:
            # Tab2 proportions
            for j in range(n_categories):
                opacities.append(1.0 if j == i else 0.15)
            
            # Tab2 counts
            for j in range(n_categories):
                opacities.append(1.0 if j == i else 0.15)
        
        dropdown_buttons.append(
            dict(
                label=category,
                method="restyle",
                args=[{"marker.opacity": opacities}]
            )
        )
    
    # Create visibility arrays based on whether tab2 exists
    if has_tab2:
        # Tab 1 Proportions
        tab1_prop_visible = [True]*n_categories + [True]*n_categories + [False]*n_categories + [False]*n_categories + [False]*n_categories
        # Tab 1 Counts
        tab1_count_visible = [True]*n_categories + [False]*n_categories + [True]*n_categories + [False]*n_categories + [False]*n_categories
        # Tab 2 Proportions
        tab2_prop_visible = [True]*n_categories + [False]*n_categories + [False]*n_categories + [True]*n_categories + [False]*n_categories
        # Tab 2 Counts
        tab2_count_visible = [True]*n_categories + [False]*n_categories + [False]*n_categories + [False]*n_categories + [True]*n_categories
    else:
        # Single tab mode
        # Proportions
        prop_visible = [True]*n_categories + [True]*n_categories + [False]*n_categories
        # Counts
        count_visible = [True]*n_categories + [False]*n_categories + [True]*n_categories
    
    # Build update menus
    updatemenus = []
    
    # Dropdown for cell type selection
    updatemenus.append(
        dict(
            buttons=dropdown_buttons,
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.35 if has_tab2 else 0.18,
            xanchor="left",
            y=1.12,
            yanchor="top",
            bgcolor="white",
            bordercolor="gray",
            borderwidth=1
        )
    )
    
    # Tab selection buttons (only if tab2 exists)
    if has_tab2:
        updatemenus.append(
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[
                            {"visible": tab1_prop_visible},
                            {"yaxis2.title": "Proportion",
                             "annotations[1].text": tab1_title}
                        ],
                        label=tab1_title,
                        method="update"
                    ),
                    dict(
                        args=[
                            {"visible": tab2_prop_visible},
                            {"yaxis2.title": "Proportion",
                             "annotations[1].text": tab2_title}
                        ],
                        label=tab2_title,
                        method="update"
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.12,
                yanchor="top"
            )
        )
    
    # Proportions/Counts toggle
    if has_tab2:
        updatemenus.append(
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[
                            {"visible": tab1_prop_visible},
                            {"yaxis2.title": "Proportion"}
                        ],
                        label="Proportions",
                        method="update"
                    ),
                    dict(
                        args=[
                            {"visible": tab1_count_visible},
                            {"yaxis2.title": "Absolute Count"}
                        ],
                        label="Absolute Counts",
                        method="update"
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.18,
                xanchor="left",
                y=1.12,
                yanchor="top"
            )
        )
    else:
        updatemenus.append(
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[
                            {"visible": prop_visible},
                            {"yaxis2.title": "Proportion"}
                        ],
                        label="Proportions",
                        method="update"
                    ),
                    dict(
                        args=[
                            {"visible": count_visible},
                            {"yaxis2.title": "Absolute Count"}
                        ],
                        label="Absolute Counts",
                        method="update"
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.12,
                yanchor="top"
            )
        )
    
    # Update layout
    fig.update_layout(
        updatemenus=updatemenus,
        barmode='stack',
        legend=dict(
            title="Cell Type",
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.02
        ),
        width=width,
        height=height,
        hovermode='closest',
        template='plotly_white'
    )
    
    # Update UMAP axes - center the plot
    umap_x_range = [umap_coords[:, 0].min() - 1, umap_coords[:, 0].max() + 1]
    umap_y_range = [umap_coords[:, 1].min() - 1, umap_coords[:, 1].max() + 1]
    
    umap_x_center = np.mean(umap_x_range)
    umap_y_center = np.mean(umap_y_range)
    umap_range = max(umap_x_range[1] - umap_x_range[0], umap_y_range[1] - umap_y_range[0]) / 2
    
    fig.update_xaxes(
        title="UMAP 1",
        showgrid=False,
        zeroline=False,
        range=[umap_x_center - umap_range, umap_x_center + umap_range],
        row=1, col=1
    )
    fig.update_yaxes(
        title="UMAP 2",
        showgrid=False,
        zeroline=False,
        scaleanchor="x",
        scaleratio=1,
        range=[umap_y_center - umap_range, umap_y_center + umap_range],
        row=1, col=1
    )
    
    # Update bar chart axes
    fig.update_xaxes(
        title="Condition",
        showgrid=False,
        row=2, col=1
    )
    fig.update_yaxes(
        title="Proportion",
        showgrid=False,
        row=2, col=1
    )
    
    # Save as self-contained interactive HTML
    fig.write_html(output_file, include_plotlyjs=True)
    
    print(f"Interactive plot saved as '{output_file}'")