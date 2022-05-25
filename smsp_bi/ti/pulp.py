import pulp

from smsp_bi import base
from smsp_bi.utils import Schedule


class TI(base.TI):
    def __init__(self, smsp, name="TI"):
        self.m = pulp.LpProblem(name, pulp.LpMinimize)
        super().__init__(smsp)

    def _update(self):
        pass

    def optimize(self):
        self.m.solve()

    def _create_x_variables(self):
        self.x_indices = [
            (j, t) for j in self.J for t in range(1, self.T - self.p[j] + 2)
        ]
        self.x_vars = pulp.LpVariable.dicts(
            "x", self.x_indices, lowBound=0, upBound=1, cat=pulp.LpInteger
        )
        self.m += pulp.lpSum(
            [self._make_cost(*ind) * self.x_vars[ind] for ind in self.x_indices]
        )

    def _add_job_completion_constraints(self):
        for job in self.J:
            self.m += (
                pulp.lpSum(self.x_vars[(j, t)] for j, t in self.x_indices if j == job)
                == 1,
                f"Completion Constraint[{job}]",
            )

    def _add_machine_capacity_constraints(self):
        for t in range(1, self.T + 1):
            self.m += (
                pulp.lpSum(
                    self.x_vars[(j, s)]
                    for j in self.J
                    for s in range(max(1, t - self.p[j] + 1), t + 1)
                    if (j, s) in self.x_indices
                )
                <= 1,
                f"Cardinality Constraint[{t}]",
            )

    def _make_start_time(self, job):
        return int(
            round(
                sum(
                    [
                        (t - 1) * self.x_vars[(j, t)].value()
                        for j, t in self.x_indices
                        if j == job
                    ]
                )
            )
        )

    def get_schedule(self):
        start_times = [self._make_start_time(j) for j in self.J]

        return Schedule(
            start_times=start_times,
            end_times=[s + p for s, p in zip(start_times, self.p)],
        )
