"""
Plot a grouped execution timeline from a TAU Chrome/Perfetto JSON trace.

The script reads TAU begin/end trace events, reconstructs execution
intervals for selected high-level CalTop timers, and displays them as a
Gantt-style horizontal timeline.

Only selected high-level CalTop timers are shown. Automatically generated
OpenMP/MKL events and lower-level Density Filter timers are excluded.
"""

import json
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns


# ============================================================
# User settings
# ============================================================

TRACE_FILE = "Threads_1/trace.json"

FONT_NAME = "Times New Roman"
FONT_SIZE = 22

DPI = 300

OUTPUT_FILE = "CalTop_trace.png"


# ============================================================
# Seaborn / Matplotlib formatting
# ============================================================

# Seaborn colorblind style
sns.set_theme(
    style="ticks",
    palette="colorblind"
)

# Times New Roman
# These settings are applied AFTER sns.set_theme() so that
# Seaborn does not overwrite them.
plt.rcParams["font.family"] = FONT_NAME
plt.rcParams["font.serif"] = [FONT_NAME]
plt.rcParams["font.size"] = FONT_SIZE

plt.rcParams["xtick.labelsize"] = FONT_SIZE
plt.rcParams["ytick.labelsize"] = FONT_SIZE

plt.rcParams["mathtext.fontset"] = "stix"


# ============================================================
# TAU timers to extract
# ============================================================

# These names must exactly match the timer names in the TAU trace.
TIMERS = [
    "preprocess: CalTop",
    "filter: CalTop",
    "Linear system: Assembly",
    "PARDISO: Factorization",
    "PARDISO: Solve",
    "PARDISO: Cleanup",
    "Stress Grad: Total",
    "Compliance Grad: Total",
    "CG Grad: Total",
    "Vol Grad: Total",
    "fileIO: CalTop",
]


# ============================================================
# Display names used in the figure
# ============================================================

DISPLAY_LABELS = {
    "preprocess: CalTop": "Preprocess()",
    "filter: CalTop": "Filter()",
    "Linear system: Assembly": "Assembly()",

    "PARDISO: Factorization": "Factorization",
    "PARDISO: Solve": "Solve",
    "PARDISO: Cleanup": "Cleanup",

    "Stress Grad: Total": "Stress()",
    "Compliance Grad: Total": "Compliance()",
    "CG Grad: Total": "Center of Gravity()",
    "Vol Grad: Total": "Volume Fraction()",

    "fileIO: CalTop": "Output()",
}


# ============================================================
# Timer groups
# ============================================================

GROUPS = {
    "Preprocessing": [
        "preprocess: CalTop",
    ],

    "Filter": [
        "filter: CalTop",
    ],

    "Linear System": [
        "Linear system: Assembly",
    ],

    "PARDISO": [
        "PARDISO: Factorization",
        "PARDISO: Solve",
        "PARDISO: Cleanup",
    ],

    "Sensitivity": [
        "Stress Grad: Total",
        "Compliance Grad: Total",
        "CG Grad: Total",
        "Vol Grad: Total",
    ],

    "File I/O": [
        "fileIO: CalTop",
    ],
}


# ============================================================
# Map each timer to its group
# ============================================================

timer_to_group = {}

for group, timers in GROUPS.items():

    for timer in timers:

        timer_to_group[timer] = group


# ============================================================
# Read TAU trace
# ============================================================

events = []

with open(TRACE_FILE, "r") as f:

    for line in f:

        line = line.strip()

        # TAU Chrome JSON may contain one event per line.
        if not line.startswith('{"name"'):
            continue

        # Remove trailing comma.
        if line.endswith(","):
            line = line[:-1]

        events.append(
            json.loads(line)
        )


# ============================================================
# Determine application start time
# ============================================================

application_start = None

for event in events:

    name = event.get(
        "name",
        ""
    ).strip()

    phase = event.get("ph")

    if (
        name == ".TAU application"
        and phase == "B"
        and str(event.get("tid", "0")) == "0"
    ):

        application_start = float(
            event["ts"]
        )

        break


if application_start is None:

    raise RuntimeError(
        "Could not find the '.TAU application' start event."
    )


# ============================================================
# Reconstruct begin/end timer intervals
# ============================================================

active_timers = {}

intervals = []


for event in events:

    name = event.get(
        "name",
        ""
    ).strip()

    # Ignore timers not included in TIMERS.
    if name not in TIMERS:
        continue

    pid = str(
        event.get("pid", "0")
    )

    tid = str(
        event.get("tid", "0")
    )

    key = (
        name,
        pid,
        tid,
    )

    timestamp = float(
        event["ts"]
    )

    phase = event.get("ph")


    # --------------------------------------------------------
    # Timer begin
    # --------------------------------------------------------

    if phase == "B":

        if key not in active_timers:

            active_timers[key] = []

        active_timers[key].append(
            timestamp
        )


    # --------------------------------------------------------
    # Timer end
    # --------------------------------------------------------

    elif phase == "E":

        if (
            key in active_timers
            and active_timers[key]
        ):

            start = active_timers[key].pop()

            start_seconds = (
                start - application_start
            ) / 1.0e6

            end_seconds = (
                timestamp - application_start
            ) / 1.0e6

            duration_seconds = (
                timestamp - start
            ) / 1.0e6

            intervals.append(
                {
                    "name": name,
                    "start": start_seconds,
                    "end": end_seconds,
                    "duration": duration_seconds,
                }
            )


# ============================================================
# Create figure
# ============================================================

fig, ax = plt.subplots()


# ============================================================
# Increase figure size
# ============================================================

F = plt.gcf()

Size = F.get_size_inches()

F.set_size_inches(
    Size[0] * 2.5,
    Size[1] * 2.0,
    forward=True
)


# ============================================================
# Seaborn colorblind palette
# ============================================================

colors = sns.color_palette(
    "colorblind",
    n_colors=len(GROUPS),
)

group_colors = {}

for i, group in enumerate(GROUPS):

    group_colors[group] = colors[i]


# ============================================================
# Plot timer intervals
# ============================================================

for i, timer in enumerate(TIMERS):

    for interval in intervals:

        if interval["name"] != timer:
            continue

        group = timer_to_group[timer]

        ax.barh(
            i,
            interval["duration"],
            left=interval["start"],
            height=0.62,
            color=group_colors[group],
            edgecolor="black",
            linewidth=0.7,
        )


# ============================================================
# Y-axis formatting
# ============================================================

ax.set_yticks(
    range(len(TIMERS))
)

ax.set_yticklabels(
    [
        DISPLAY_LABELS[timer]
        for timer in TIMERS
    ],
    fontsize=FONT_SIZE,
    fontname=FONT_NAME,
)

# Remove Y-axis title
ax.set_ylabel("")

ax.invert_yaxis()


# ============================================================
# X-axis formatting
# ============================================================

ax.set_xlabel(
    "Execution Time (s)",
    fontsize=FONT_SIZE,
    fontname=FONT_NAME,
)

ax.tick_params(
    axis="x",
    which="major",
    labelsize=FONT_SIZE,
)

ax.tick_params(
    axis="x",
    which="minor",
    length=4,
)


# ============================================================
# Minor ticks
# ============================================================

# One minor tick halfway between adjacent major ticks.
ax.xaxis.set_minor_locator(
    AutoMinorLocator(2)
)

# Disable minor ticks on categorical Y-axis.
ax.tick_params(
    axis="y",
    which="minor",
    left=False,
)


# ============================================================
# Explicitly set X-axis tick-label font
# ============================================================

# Major X-axis tick labels
for label in ax.get_xticklabels(
    which="major"
):

    label.set_fontname(
        FONT_NAME
    )

    label.set_fontsize(
        FONT_SIZE
    )


# Minor X-axis tick labels, if displayed
for label in ax.get_xticklabels(
    which="minor"
):

    label.set_fontname(
        FONT_NAME
    )

    label.set_fontsize(
        FONT_SIZE
    )


# ============================================================
# Explicitly set Y-axis tick-label font
# ============================================================

for label in ax.get_yticklabels():

    label.set_fontname(
        FONT_NAME
    )

    label.set_fontsize(
        FONT_SIZE
    )


# ============================================================
# Major grid
# ============================================================

ax.grid(
    axis="x",
    which="major",
    linestyle="--",
    linewidth=0.8,
    alpha=1.0,
)


# ============================================================
# Minor grid
# ============================================================

ax.grid(
    axis="x",
    which="minor",
    linestyle="--",
    linewidth=0.8,
    alpha=1.0,
)

ax.set_axisbelow(True)


# ============================================================
# Draw separators between timer groups
# ============================================================

previous_group = None

for i, timer in enumerate(TIMERS):

    current_group = timer_to_group[timer]

    if (
        previous_group is not None
        and current_group != previous_group
    ):

        ax.axhline(
            i - 0.5,
            linewidth=0.5,
            alpha=0.30,
        )

    previous_group = current_group


# ============================================================
# X-axis limits
# ============================================================

if intervals:

    maximum_time = max(
        interval["end"]
        for interval in intervals
    )

    ax.set_xlim(
        0.0,
        maximum_time * 1.02,
    )


# ============================================================
# Legend
# ============================================================

legend_handles = []

for group in GROUPS:

    legend_handles.append(
        Patch(
            facecolor=group_colors[group],
            edgecolor="black",
            label=group,
        )
    )


legend = ax.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.16),
    ncol=3,
    frameon=False,
    fontsize=FONT_SIZE,
)


# Explicitly set legend font
for text in legend.get_texts():

    text.set_fontname(
        FONT_NAME
    )

    text.set_fontsize(
        FONT_SIZE
    )


# ============================================================
# Final formatting
# ============================================================

plt.tight_layout()


# ============================================================
# Save figure
# ============================================================

plt.savefig(
    OUTPUT_FILE,
    dpi=DPI,
    bbox_inches="tight",
)


# ============================================================
# Display figure
# ============================================================

plt.show()