from abc import ABC, abstractmethod


class TI(ABC):
    def __init__(self, smsp):
        # smsp is SMSP object
        self.smsp = smsp.copy()
        self.p = smsp._processing_times
        self.d = smsp._due_dates
        self.c = smsp._cost
        self.T = sum(self.p)
        self.J = range(len(self.p))
        self._create_model()

    def _make_cost(self, j, t):
        return self.c[j] * max(0, t - 1 - self.d[j])

    @abstractmethod
    def _update(self):
        pass

    @abstractmethod
    def _create_x_variables(self):
        pass

    @abstractmethod
    def _add_job_completion_constraints(self):
        pass

    @abstractmethod
    def _add_machine_capacity_constraints(self):
        pass

    def _create_model(self):
        self._create_x_variables()
        self._update()
        self._add_job_completion_constraints()
        self._add_machine_capacity_constraints()

    @abstractmethod
    def optimize(self):
        pass

    @abstractmethod
    def get_schedule(self):
        pass
