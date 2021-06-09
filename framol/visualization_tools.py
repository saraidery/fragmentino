class VisualizationTool:
    def update_figure_layout(self, fig):

        fig.update_layout(
            showlegend=False,
            scene=dict(
                xaxis_title="",
                yaxis_title="",
                zaxis_title="",
                xaxis=dict(
                    showbackground=False,
                    tickvals=[],
                ),
                yaxis=dict(
                    showbackground=False,
                    tickvals=[],
                ),
                zaxis=dict(
                    showbackground=False,
                    tickvals=[],
                ),
            ),
            margin=dict(l=0, r=0, t=0, b=0),
        )

    def show_figure(self, fig):
        fig.show()  # pragma: no cover
