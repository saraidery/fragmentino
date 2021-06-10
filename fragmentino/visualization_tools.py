#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import plotly.graph_objects as go


class Figure:
    def __init__(self, data):
        self.fig = go.Figure(data)

    def update_figure_layout(self):
        self.fig.update_layout(
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

    def show_figure(self, **kwargs):
        self.fig.show(**kwargs)  # pragma: no cover
