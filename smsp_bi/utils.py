import numpy as np


def irange(a, b):
    return range(a, b + 1)


class Schedule:
    def __init__(self, start_times, end_times):
        self.start_times = np.array(start_times)
        self.end_times = np.array(end_times)

    def validate(self):
        start_times = sorted(self.start_times)
        end_times = self.end_times[self.start_times.argsort()]
        assert len((start_times[1:] - end_times[:-1] < 0).nonzero()[0]) == 0


class SMSP:
    def __init__(self, processing_times, due_dates=None, cost=None):

        if due_dates is None:
            due_dates = [0] * len(processing_times)
        if cost is None:
            cost = [1] * len(processing_times)

        self._processing_times = np.array(processing_times)
        self._due_dates = np.array(due_dates)
        self._cost = np.array(cost)

    def get_schedule_from_sequence(self, sequence):
        start_times = np.cumsum([p[j] for j in sequence])[np.argsort(sequence)]
        end_times = start_times + self._processing_times
        return Schedule(start_times, end_times)

    def get_objective_from_schedule(self, schedule):
        due_delta = schedule.start_times - self._due_dates
        return np.where(due_delta > 0, due_delta, 0) @ self._cost

    def get_objective_from_sequence(self, sequence):
        return self.get_objective_from_schedule(
            self.get_schedule_from_sequence(sequence)
        )

    def print_specs(self):
        print(f"Number of jobs: {len(self._processing_times)}")
        print(f"Processing times: {list(self._processing_times)}")
        print(f"Due dates: {list(self._due_dates)}")
        print(f"Tardy Cost: {list(self._cost)}")

    def copy(self):
        return SMSP(
            self._processing_times,
            self._due_dates,
            self._processing_times,
        )


def get_example_problem():
    return SMSP(
        processing_times=[10, 11, 14, 11, 12, 10, 10, 12, 10, 14, 13, 14, 14],
        due_dates=[24, 22, 37, 20, 37, 31, 35, 27, 28, 39, 39, 29, 22],
        cost=[15, 18, 18, 12, 20, 11, 7, 11, 21, 12, 7, 22, 7],
    )
