import gurobipy as gp
from gurobipy import GRB

from smsp_bi import base
from smsp_bi.utils import Schedule


class TI(base.TI):
    def __init__(self, smsp, name="TI"):
        self.m = gp.Model(name)
        self.m.setAttr("ModelSense", GRB.MINIMIZE)
        super().__init__(smsp)

    def _update(self):
        self.m.update()

    def optimize(self):
        self.m.optimize()

    def _create_x_variables(self):
        self.x_indices, cost = gp.multidict(
            {
                (j, t): self._make_cost(j, t)
                for j in self.J
                for t in range(1, self.T - self.p[j] + 2)
            }
        )
        self.x_vars = self.m.addVars(
            self.x_indices, obj=cost, vtype=GRB.BINARY, name="x"
        )

    def _add_job_completion_constraints(self):
        self.m.addConstrs(self.x_vars.sum(j, "*") == 1 for j in self.J)

    def _add_machine_capacity_constraints(self):
        self.m.addConstrs(
            gp.quicksum(
                self.x_vars[(j, s)]
                for j in self.J
                for s in range(max(1, t - self.p[j] + 1), t + 1)
                if (j, s) in self.x_indices
            )
            <= 1
            for t in range(1, self.T + 1)
        )

    def _make_start_time(self, job):
        return int(
            round(
                gp.quicksum(
                    (t - 1) * self.x_vars[(j, t)] for j, t in self.x_indices if j == job
                ).getValue()
            )
        )

    def get_schedule(self):
        start_times = [self._make_start_time(j) for j in self.J]

        return Schedule(
            start_times=start_times,
            end_times=[s + p for s, p in zip(start_times, self.p)],
        )
