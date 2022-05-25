import gurobipy as gp
from gurobipy import GRB

from smsp_bi import base
from smsp_bi.utils import Schedule, irange


class BI(base.BI):
    def __init__(self, smsp, name="BI"):
        self.m = gp.Model(name)
        self.m.setAttr("ModelSense", GRB.MINIMIZE)
        super().__init__(smsp)

    def _update(self):
        self.m.update()

    def optimize(self):
        self.m.optimize()

    def _create_z_u_variables(self):
        self.z_vars = self.m.addVars(self.z_indices, obj=0, vtype=GRB.BINARY, name="z")
        self.u_vars = self.m.addVars(
            self.z_indices, obj=0, vtype=GRB.CONTINUOUS, name="u"
        )

    def _add_job_completion_constraints(self):
        self.m.addConstrs(self.z_vars.sum(j, "*", "*") == 1 for j in self.J)

    def _add_machine_capacity_constraints_1(self):
        self.m.addConstrs(
            gp.quicksum(
                self.z_vars[(j, a, k)]
                for j in self.J
                for k in self.K[j]
                for a in irange(max(1, b - self.P[j] - k + 2), b)
                if (j, a, k) in self.z_indices
            )
            <= 1
            for b in irange(1, self.B)
        )

    def _add_machine_capacity_constraints_2(self):
        self.m.addConstrs(
            self.u_vars.sum("*", b, "*")
            - gp.quicksum(
                self.u_vars[(j, b - self.P[j] - k + 1, k)]
                for j in self.J
                for k in self.K[j]
                if (j, b - self.P[j] - k + 1, k) in self.z_indices
            )
            + gp.quicksum(
                (2 - k - self.pi[j]) * self.z_vars[(j, b - self.P[j] - k + 1, k)]
                for j in self.J
                for k in self.K[j]
                if (j, b - self.P[j] - k + 1, k) in self.z_indices
            )
            + gp.quicksum(
                self.z_vars[(j, a, k)]
                for j in self.J
                for k in self.K[j]
                for a in irange(max(1, b - self.P[j] - k + 2), b - 1)
                if (j, a, k) in self.z_indices
            )
            <= 1
            for b in irange(1, self.B)
        )

    def _create_u_lower_bound_constraints(self):
        self.m.addConstrs(
            ((1 - k) * (1 - self.pi[j]) + 1 / self.Delta) * self.z_vars[(j, b, k)]
            - self.u_vars[(j, b, k)]
            <= 0
            for j, b, k in self.z_indices
        )

    def _create_u_upper_bound_constraints(self):
        self.m.addConstrs(
            self.u_vars[(j, b, k)] - (1 - k * self.pi[j]) * self.z_vars[(j, b, k)] <= 0
            for j, b, k in self.z_indices
        )

    def _create_T_variables(self):
        T_indices, cost = gp.multidict(
            {
                (j, b, k): self.c[j] * self.Delta
                for j in self.J
                for k in self.K[j]
                for b in irange(self.D[j], self.B - self.P[j] + 1)
            }
        )
        self.T_vars = self.m.addVars(
            T_indices, obj=cost, vtype=GRB.CONTINUOUS, name="T"
        )

    def _create_T_k0_lower_bound_constraints(self):
        self.m.addConstrs(
            self.delta[j] * self.z_vars[(j, self.D[j], 0)]
            - self.u_vars[(j, self.D[j], 0)]
            - self.T_vars[(j, self.D[j], 0)]
            <= 0
            for j in self.J
            if 1 - self.pi[j] < self.delta[j]
        )

    def _create_T_k0_upper_bound_constraints(self):
        self.m.addConstrs(
            self.T_vars[(j, self.D[j], 0)]
            - (self.delta[j] - 1 + self.pi[j]) * self.z_vars[(j, self.D[j], 0)]
            <= 0
            for j in self.J
            if 1 - self.pi[j] < self.delta[j]
        )

    def _create_T_k1_equality_constraints(self):
        self.m.addConstrs(
            self.T_vars[(j, self.D[j], 1)]
            - self.delta[j] * self.z_vars[(j, self.D[j], 1)]
            + self.u_vars[(j, self.D[j], 1)]
            == 0
            for j in self.J
            if 1 in self.K[j] and 1 - self.pi[j] < self.delta[j]
        )

    def _create_T_k1_lower_bound_constraints(self):
        self.m.addConstrs(
            self.delta[j] * self.z_vars[(j, self.D[j], 1)]
            - self.u_vars[(j, self.D[j], 1)]
            - self.T_vars[(j, self.D[j], 1)]
            <= 0
            for j in self.J
            if 1 in self.K[j] and self.delta[j] <= 1 - self.pi[j]
        )

    def _create_T_k1_upper_bound_constraints(self):
        self.m.addConstrs(
            self.T_vars[(j, self.D[j], 1)]
            - self.delta[j] * self.z_vars[(j, self.D[j], 1)]
            <= 0
            for j in self.J
            if 1 in self.K[j] and self.delta[j] <= 1 - self.pi[j]
        )

    def _create_T_equality_constraints(self):
        self.m.addConstrs(
            self.T_vars[(j, b, k)]
            - (b - self.D[j] + self.delta[j]) * self.z_vars[(j, b, k)]
            + self.u_vars[(j, b, k)]
            == 0
            for j in self.J
            for k in self.K[j]
            for b in irange(self.D[j] + 1, self.B)
            if (j, b, k) in self.z_indices
        )

    def _make_start_time(self, job):
        return int(
            round(
                gp.quicksum(
                    self.Delta * (b * self.z_vars[(j, b, k)] - self.u_vars[(j, b, k)])
                    for j, b, k in self.z_indices
                    if j == job
                ).getValue()
            )
        )

    def get_schedule(self):
        start_times = [self._make_start_time(j) for j in self.J]

        return Schedule(
            start_times=start_times,
            end_times=[s + p for s, p in zip(start_times, self.p)],
        )


class BI_slim(BI):
    def _setup(self):
        super()._setup()
        self._setup_slim()

    def _add_machine_capacity_constraints_2(self):
        self.m.addConstrs(
            self.Delta * self.u_vars.sum("*", b, "*")
            - gp.quicksum(
                self.Delta * self.u_vars[(j, b - self.P[j] - k + 1, k)]
                for j in self.J
                for k in self.K[j]
                if (j, b - self.P[j] - k + 1, k) in self.z_indices
            )
            + gp.quicksum(
                (2 * self.Delta - self.Delta * k - self.Delta_pi[j])
                * self.z_vars[(j, b - self.P[j] - k + 1, k)]
                for j in self.J
                for k in self.K[j]
                if (j, b - self.P[j] - k + 1, k) in self.z_indices
            )
            + gp.quicksum(
                self.Delta * self.z_vars[(j, a, k)]
                for j in self.J
                for k in self.K[j]
                for a in irange(max(1, b - self.P[j] - k + 2), b - 1)
                if (j, a, k) in self.z_indices
            )
            <= self.Delta
            for b in range(1, self.B + 1)
        )

    def _create_u_lower_bound_constraints(self):
        self.m.addConstrs(
            ((1 - k) * (self.Delta - self.Delta_pi[j]) + 1) * self.z_vars[(j, b, k)]
            - self.Delta * self.u_vars[(j, b, k)]
            <= 0
            for j, b, k in self.z_indices
        )

    def _create_u_upper_bound_constraints(self):
        self.m.addConstrs(
            self.Delta * self.u_vars[(j, b, k)]
            - (self.Delta - k * self.Delta_pi[j]) * self.z_vars[(j, b, k)]
            <= 0
            for j, b, k in self.z_indices
        )

    def _create_T_k0_lower_bound_constraints(self):
        self.m.addConstrs(
            self.Delta_delta[j] * self.z_vars[(j, self.D[j], 0)]
            - self.Delta * self.u_vars[(j, self.D[j], 0)]
            - self.Delta * self.T_vars[(j, self.D[j], 0)]
            <= 0
            for j in self.J
            if self.Delta - self.Delta_pi[j] < self.Delta_delta[j]
        )

    def _create_T_k0_upper_bound_constraints(self):
        self.m.addConstrs(
            self.Delta * self.T_vars[(j, self.D[j], 0)]
            - (self.Delta_delta[j] - self.Delta + self.Delta_pi[j])
            * self.z_vars[(j, self.D[j], 0)]
            <= 0
            for j in self.J
            if self.Delta - self.Delta_pi[j] < self.Delta_delta[j]
        )

    def _create_T_k1_equality_constraints(self):
        self.m.addConstrs(
            self.Delta * self.T_vars[(j, self.D[j], 1)]
            - self.Delta_delta[j] * self.z_vars[(j, self.D[j], 1)]
            + self.Delta * self.u_vars[(j, self.D[j], 1)]
            == 0
            for j in self.J
            if 1 in self.K[j] and self.Delta - self.Delta_pi[j] < self.Delta_delta[j]
        )

    def _create_T_k1_lower_bound_constraints(self):
        self.m.addConstrs(
            self.Delta_delta[j] * self.z_vars[(j, self.D[j], 1)]
            - self.Delta * self.u_vars[(j, self.D[j], 1)]
            - self.Delta * self.T_vars[(j, self.D[j], 1)]
            <= 0
            for j in self.J
            if 1 in self.K[j] and self.Delta_delta[j] <= self.Delta - self.Delta_pi[j]
        )

    def _create_T_k1_upper_bound_constraints(self):
        self.m.addConstrs(
            self.Delta * self.T_vars[(j, self.D[j], 1)]
            - self.Delta_delta[j] * self.z_vars[(j, self.D[j], 1)]
            <= 0
            for j in self.J
            if 1 in self.K[j] and self.Delta_delta[j] <= self.Delta - self.Delta_pi[j]
        )

    def _create_T_equality_constraints(self):
        self.m.addConstrs(
            self.Delta * self.T_vars[(j, b, k)]
            - (self.Delta * b - self.Delta * self.D[j] + self.Delta_delta[j])
            * self.z_vars[(j, b, k)]
            + self.Delta * self.u_vars[(j, b, k)]
            == 0
            for j in self.J
            for k in self.K[j]
            for b in irange(self.D[j] + 1, self.B)
            if (j, b, k) in self.z_indices
        )
