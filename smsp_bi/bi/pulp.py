import pulp

from smsp_bi import base
from smsp_bi.utils import Schedule, irange


class BI_2(base.BI_2):
    def __init__(self, smsp, name="BI_2"):
        self.m = pulp.LpProblem(name, pulp.LpMinimize)
        super().__init__(smsp)

    def _update(self):
        pass

    def optimize(self):
        self.m.solve()

    def _create_z_u_variables(self):
        self.z_vars = pulp.LpVariable.dicts(
            "z", self.z_indices, lowBound=0, upBound=1, cat=pulp.LpInteger
        )
        self.u_vars = pulp.LpVariable.dicts(
            "u", self.z_indices, lowBound=0, cat=pulp.LpContinuous
        )

    def _add_job_completion_constraints(self):
        for job in self.J:
            self.m += (
                pulp.lpSum(
                    self.z_vars[(j, b, k)] for j, b, k in self.z_indices if j == job
                )
                == 1,
                f"Completion Constraint[{job}]",
            )

    def _add_machine_capacity_constraints_1(self):
        for b in irange(1, self.B):
            self.m += (
                pulp.lpSum(
                    self.z_vars[(j, a, k)]
                    for j in self.J
                    for k in self.K[j]
                    for a in irange(max(1, b - self.P[j] - k + 2), b)
                    if (j, a, k) in self.z_indices
                )
                <= 1,
                f"Completion Constraints1 [{b}]",
            )

    def _add_machine_capacity_constraints_2(self):
        for b in irange(1, self.B):
            self.m += (
                pulp.lpSum(
                    self.u_vars[(j, a, k)] for (j, a, k) in self.z_indices if a == b
                )
                - pulp.lpSum(
                    self.u_vars[(j, b - self.P[j] - k + 1, k)]
                    for j in self.J
                    for k in self.K[j]
                    if (j, b - self.P[j] - k + 1, k) in self.z_indices
                )
                + pulp.lpSum(
                    (2 - k - self.pi[j]) * self.z_vars[(j, b - self.P[j] - k + 1, k)]
                    for j in self.J
                    for k in self.K[j]
                    if (j, b - self.P[j] - k + 1, k) in self.z_indices
                )
                + pulp.lpSum(
                    self.z_vars[(j, a, k)]
                    for j in self.J
                    for k in self.K[j]
                    for a in irange(max(1, b - self.P[j] - k + 2), b - 1)
                    if (j, a, k) in self.z_indices
                )
                <= 1,
                f"Completion Constraints2 [{b}]",
            )

    def _create_u_lower_bound_constraints(self):
        for j, b, k in self.z_indices:
            self.m += (
                ((1 - k) * (1 - self.pi[j]) + 1 / self.Delta) * self.z_vars[(j, b, k)]
                - self.u_vars[(j, b, k)]
                <= 0,
                f"u lower bound constraints [{j, b, k}]",
            )

    def _create_u_upper_bound_constraints(self):
        for j, b, k in self.z_indices:
            self.m += (
                self.u_vars[(j, b, k)] - (1 - k * self.pi[j]) * self.z_vars[(j, b, k)]
                <= 0,
                f"u upper bound constraints [{j, b, k}]",
            )

    def _create_T_variables(self):
        self.T_indices = [
            (j, b, k)
            for j in self.J
            for k in self.K[j]
            for b in irange(self.D[j], self.B - self.P[j] + 1)
        ]
        self.T_vars = pulp.LpVariable.dicts(
            "T", self.T_indices, lowBound=0, cat=pulp.LpContinuous
        )

        self.m += pulp.lpSum(
            self.c[j] * self.Delta * self.T_vars[(j, b, k)]
            for j, b, k in self.T_indices
        )

    def _create_T_k0_lower_bound_constraints(self):
        for j in self.J:
            if 1 - self.pi[j] < self.delta[j]:
                self.m += (
                    self.delta[j] * self.z_vars[(j, self.D[j], 0)]
                    - self.u_vars[(j, self.D[j], 0)]
                    - self.T_vars[(j, self.D[j], 0)]
                    <= 0,
                    f"T_k0 lower bound constraints [{j}]",
                )

    def _create_T_k0_upper_bound_constraints(self):
        for j in self.J:
            if 1 - self.pi[j] < self.delta[j]:
                self.m += (
                    self.T_vars[(j, self.D[j], 0)]
                    - (self.delta[j] - 1 + self.pi[j]) * self.z_vars[(j, self.D[j], 0)]
                    <= 0,
                    f"T_k0 upper bound constraints [{j}]",
                )

    def _create_T_k1_equality_constraints(self):
        for j in self.J:
            if 1 in self.K[j] and 1 - self.pi[j] < self.delta[j]:
                self.m += (
                    self.T_vars[(j, self.D[j], 1)]
                    - self.delta[j] * self.z_vars[(j, self.D[j], 1)]
                    + self.u_vars[(j, self.D[j], 1)]
                    == 0,
                    f"T_k1 equality constraints [{j}]",
                )

    def _create_T_k1_lower_bound_constraints(self):
        for j in self.J:
            if 1 in self.K[j] and self.delta[j] <= 1 - self.pi[j]:
                self.m += (
                    self.delta[j] * self.z_vars[(j, self.D[j], 1)]
                    - self.u_vars[(j, self.D[j], 1)]
                    - self.T_vars[(j, self.D[j], 1)]
                    <= 0,
                    f"T_k1 lower bound constraints [{j}]",
                )

    def _create_T_k1_upper_bound_constraints(self):
        for j in self.J:
            if 1 in self.K[j] and self.delta[j] <= 1 - self.pi[j]:
                self.m += (
                    self.T_vars[(j, self.D[j], 1)]
                    - self.delta[j] * self.z_vars[(j, self.D[j], 1)]
                    <= 0,
                    f"T_k1 upper bound constraints [{j}]",
                )

    def _create_T_equality_constraints(self):
        for j, b, k in self.T_indices:
            if b >= self.D[j] + 1:
                self.m += (
                    self.T_vars[(j, b, k)]
                    - (b - self.D[j] + self.delta[j]) * self.z_vars[(j, b, k)]
                    + self.u_vars[(j, b, k)]
                    == 0,
                    f"T_equality constraints [{j,b,k}]",
                )

    def _make_start_time(self, job):
        return int(
            round(
                sum(
                    [
                        self.Delta
                        * (
                            b * self.z_vars[(j, b, k)].value()
                            - self.u_vars[(j, b, k)].value()
                        )
                        for j, b, k in self.z_indices
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


class BI_slim(BI):
    def _setup(self):
        super()._setup()
        self._setup_slim()

    def _add_machine_capacity_constraints_2(self):
        for b in irange(1, self.B):
            self.m += (
                pulp.lpSum(
                    self.Delta * self.u_vars[(j, a, k)]
                    for (j, a, k) in self.z_indices
                    if a == b
                )
                - pulp.lpSum(
                    self.Delta * self.u_vars[(j, b - self.P[j] - k + 1, k)]
                    for j in self.J
                    for k in self.K[j]
                    if (j, b - self.P[j] - k + 1, k) in self.z_indices
                )
                + pulp.lpSum(
                    (2 * self.Delta - self.Delta * k - self.Delta_pi[j])
                    * self.z_vars[(j, b - self.P[j] - k + 1, k)]
                    for j in self.J
                    for k in self.K[j]
                    if (j, b - self.P[j] - k + 1, k) in self.z_indices
                )
                + pulp.lpSum(
                    self.Delta * self.z_vars[(j, a, k)]
                    for j in self.J
                    for k in self.K[j]
                    for a in irange(max(1, b - self.P[j] - k + 2), b - 1)
                    if (j, a, k) in self.z_indices
                )
                <= self.Delta,
                f"Completion Constraints2 [{b}]",
            )

    def _create_u_lower_bound_constraints(self):
        for j, b, k in self.z_indices:
            self.m += (
                ((1 - k) * (self.Delta - self.Delta_pi[j]) + 1) * self.z_vars[(j, b, k)]
                - self.Delta * self.u_vars[(j, b, k)]
                <= 0,
                f"u lower bound constraints [{j, b, k}]",
            )

    def _create_u_upper_bound_constraints(self):
        for j, b, k in self.z_indices:
            self.m += (
                self.Delta * self.u_vars[(j, b, k)]
                - (self.Delta - k * self.Delta_pi[j]) * self.z_vars[(j, b, k)]
                <= 0,
                f"u upper bound constraints [{j, b, k}]",
            )

    def _create_T_k0_lower_bound_constraints(self):
        for j in self.J:
            if self.Delta - self.Delta_pi[j] < self.Delta_delta[j]:
                self.m += (
                    self.Delta_delta[j] * self.z_vars[(j, self.D[j], 0)]
                    - self.Delta * self.u_vars[(j, self.D[j], 0)]
                    - self.Delta * self.T_vars[(j, self.D[j], 0)]
                    <= 0,
                    f"T_k0 lower bound constraints [{j}]",
                )

    def _create_T_k0_upper_bound_constraints(self):
        for j in self.J:
            if self.Delta - self.Delta_pi[j] < self.Delta_delta[j]:
                self.m += (
                    self.Delta * self.T_vars[(j, self.D[j], 0)]
                    - (self.Delta_delta[j] - self.Delta + self.Delta_pi[j])
                    * self.z_vars[(j, self.D[j], 0)]
                    <= 0,
                    f"T_k0 upper bound constraints [{j}]",
                )

    def _create_T_k1_equality_constraints(self):
        for j in self.J:
            if 1 in self.K[j] and self.Delta - self.Delta_pi[j] < self.Delta_delta[j]:
                self.m += (
                    self.Delta * self.T_vars[(j, self.D[j], 1)]
                    - self.Delta_delta[j] * self.z_vars[(j, self.D[j], 1)]
                    + self.Delta * self.u_vars[(j, self.D[j], 1)]
                    == 0,
                    f"T_k1 equality constraints [{j}]",
                )

    def _create_T_k1_lower_bound_constraints(self):
        for j in self.J:
            if 1 in self.K[j] and self.Delta_delta[j] <= self.Delta - self.Delta_pi[j]:
                self.m += (
                    self.Delta_delta[j] * self.z_vars[(j, self.D[j], 1)]
                    - self.Delta * self.u_vars[(j, self.D[j], 1)]
                    - self.Delta * self.T_vars[(j, self.D[j], 1)]
                    <= 0,
                    f"T_k1 lower bound constraints [{j}]",
                )

    def _create_T_k1_upper_bound_constraints(self):
        for j in self.J:
            if 1 in self.K[j] and self.Delta_delta[j] <= self.Delta - self.Delta_pi[j]:
                self.m += (
                    self.Delta * self.T_vars[(j, self.D[j], 1)]
                    - self.Delta_delta[j] * self.z_vars[(j, self.D[j], 1)]
                    <= 0,
                    f"T_k1 upper bound constraints [{j}]",
                )

    def _create_T_equality_constraints(self):
        for j, b, k in self.T_indices:
            if b >= self.D[j] + 1:
                self.m += (
                    self.Delta * self.T_vars[(j, b, k)]
                    - (self.Delta * b - self.Delta * self.D[j] + self.Delta_delta[j])
                    * self.z_vars[(j, b, k)]
                    + self.Delta * self.u_vars[(j, b, k)]
                    == 0,
                    f"T_equality constraints [{j,b,k}]",
                )
