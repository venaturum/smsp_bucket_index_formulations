import math
from abc import ABC, abstractmethod

from smsp_bi.utils import irange


class BI_2(ABC):
    def __init__(self, smsp):
        # smsp is SMSP object
        self.smsp = smsp.copy()
        self.p = smsp._processing_times
        self.d = smsp._due_dates
        self.c = smsp._cost

        self._setup()
        self._create_model()

    def _setup(self):
        T = sum(self.p)
        self.J = range(len(self.p))
        self.Delta = min(self.p)
        self.B = int(math.ceil(T / self.Delta))
        self.P = [p_ // self.Delta + 1 for p_ in self.p]
        self.pi = [P_ - p_ / self.Delta for P_, p_ in zip(self.P, self.p)]
        self.K = [(0, 1) if p_ % self.Delta else (0,) for p_ in self.p]
        self.D = [d_ // self.Delta + 1 for d_ in self.d]
        self.delta = [D_ - d_ / self.Delta for D_, d_ in zip(self.D, self.d)]
        self.z_indices = [
            (j, b, k)
            for j in self.J
            for k in self.K[j]
            for b in irange(1, self.B - self.P[j] + 1)
        ]

    def _setup_slim(self):
        self.Delta_pi = [self.Delta * P_ - p_ for P_, p_ in zip(self.P, self.p)]
        self.Delta_delta = [self.Delta * D_ - d_ for D_, d_ in zip(self.D, self.d)]

    @abstractmethod
    def _update(self):
        pass

    @abstractmethod
    def _create_z_u_variables(self):
        pass

    @abstractmethod
    def _add_job_completion_constraints(self):
        pass

    @abstractmethod
    def _add_machine_capacity_constraints_1(self):
        pass

    @abstractmethod
    def _add_machine_capacity_constraints_2(self):
        pass

    @abstractmethod
    def _create_u_lower_bound_constraints(self):
        pass

    @abstractmethod
    def _create_u_upper_bound_constraints(self):
        pass

    def _create_base_model(self):
        self._create_z_u_variables()
        self._update()
        self._add_job_completion_constraints()
        self._add_machine_capacity_constraints_1()
        self._add_machine_capacity_constraints_2()
        self._create_u_lower_bound_constraints()
        self._create_u_upper_bound_constraints()

    @abstractmethod
    def _create_T_variables(self):
        pass

    @abstractmethod
    def _create_T_k0_lower_bound_constraints(self):
        pass

    @abstractmethod
    def _create_T_k0_upper_bound_constraints(self):
        pass

    @abstractmethod
    def _create_T_k1_equality_constraints(self):
        pass

    @abstractmethod
    def _create_T_k1_lower_bound_constraints(self):
        pass

    @abstractmethod
    def _create_T_k1_upper_bound_constraints(self):
        pass

    @abstractmethod
    def _create_T_equality_constraints(self):
        pass

    def _create_tardy_vars_constraints(self):
        self._create_T_variables()
        self._update()
        self._create_T_k0_lower_bound_constraints()
        self._create_T_k0_upper_bound_constraints()
        self._create_T_k1_equality_constraints()
        self._create_T_k1_lower_bound_constraints()
        self._create_T_k1_upper_bound_constraints()
        self._create_T_equality_constraints()

    def _create_model(self):
        self._create_base_model()
        self._create_tardy_vars_constraints()
        self._update()

    @abstractmethod
    def optimize(self):
        pass

    @abstractmethod
    def get_schedule(self):
        pass
