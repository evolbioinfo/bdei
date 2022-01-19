import numpy as np
from abc import abstractmethod

EXPOSED = 'e'
INFECTED = 'i'


SAMPLING = 'sampling'
TRANSMISSION = 'transmission'
TRANSITION = 'transition'


class Model(object):
    def __init__(self, minus_avg_sigma=True, ps=None, *args, **kwargs):
        super(Model, self).__init__()
        self.__states = self._get_states()
        self.__ps = np.array(ps) if ps is not None else np.ones(len(self.states), dtype=float)
        self.__rates = np.zeros(shape=(4, len(self.states)), dtype=np.float)
        self.__sigmas = None
        self.__avg_sigma = None
        self.__minus_avg = minus_avg_sigma

    @property
    def ps(self):
        return self.__ps

    @property
    def states(self):
        return self.__states

    @property
    def rates(self):
        """
        Get rate array with states as columns,
            and transition, transmission, sampling rates and equilibrium frequencies as rows.

        :return rate array
        :rtype np.array
        """
        return self.__rates

    @abstractmethod
    def num_params(self, type=None):
        return 0

    @abstractmethod
    def _get_states(self):
        pass

    @abstractmethod
    def params2rates(self, ps, sampled_pis=None, **kwargs):
        """
        Converts parameters into a rate array.
        """
        pass

    @abstractmethod
    def get_name(self):
        pass


class State:
    def __init__(self, name, index, next_state=None, recipient=None):
        self._name = name
        self._i = index
        self._next = next_state
        self._recipient = recipient

    @property
    def index(self):
        return self._i

    @property
    def name(self):
        return self._name

    @property
    def next_state(self):
        return self._next

    @next_state.setter
    def next_state(self, next_state):
        self._next = next_state

    @property
    def recipient(self):
        return self._recipient

    @recipient.setter
    def recipient(self, recipient):
        self._recipient = recipient

    def is_symmetric(self):
        return self == self.recipient

    def __str__(self):
        return self._name

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        """Override the default equals behavior"""
        if isinstance(other, self.__class__):
            return self._name == other._name
        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

    def __hash__(self):
        """Override the default hash behavior (that returns the id or the object)"""
        return hash(self._name)


class BirthDeathExposedModel(Model):

    def __init__(self, p=0.5, *args, **kwargs):
        Model.__init__(self, ps=[0, p], minus_avg_sigma=False, *args, **kwargs)
        self.ues, self.uis = np.ones(1000, dtype=np.float128), np.ones(1000, dtype=np.float128)
        self.dt = None

    def num_params(self, type=None):
        if SAMPLING == type:
            return 1
        if TRANSMISSION == type:
            return 1
        if TRANSITION == type:
            return 1
        return 3

    def clone(self):
        model = BirthDeathExposedModel(p=self.ps[1])
        model.params2rates([self.rates[0, 0], self.rates[1, 1], self.rates[2, 1]])
        return model

    def _get_states(self):
        infected_state = State(name=INFECTED, index=1)
        exposed_state = State(name=EXPOSED, index=0, next_state=infected_state)

        # we'll put exposed as a recipient in order not to break the simulator,
        # but as the transmission rate from this state is 0, it will be ok anyway
        exposed_state.recipient = exposed_state
        infected_state.recipient = exposed_state

        return np.array([exposed_state, infected_state])

    def params2rates(self, ps, **kwargs):
        """
        Converts parameters into a rate array.

        :param ps: parameters, in the following order:
            transition E->I, transmission from I to E, sampling of I
        """
        mu, lambda_i, psi_i = ps
        self.rates[0, 0] = mu
        self.rates[1, 1] = lambda_i
        self.rates[2, 1] = psi_i

    def check_rates(self):
        return np.all(self.rates >= 0)

    def map_state(self, state):
        state = str(state)
        if state.startswith('e'):
            return self.states[0]
        if state.startswith('i'):
            return self.states[1]
        return self.states

    def get_name(self):
        return 'BDEI'
