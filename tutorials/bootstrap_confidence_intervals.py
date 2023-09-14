
# %%

from datetime import datetime, timedelta
import pickle
import sys
import numpy as np
module_path = os.path.abspath(os.path.join('../src/ocloc'))
if module_path not in sys.path:
    sys.path.append(module_path)
    from ocloc import ClockDrift, read_correlation_file, trim_correlation_trace, correlations_with_parameters
    from ocloc import ProcessingParameters, ClockDrift, suppress_stdout

    import bootstrap_module
# %%


def bootstrap_cd(clock_drift):
    """
    Bootstrap the t_app values and return a new ClockDrift object with 
    the bootstrapped values.
    """
    bootstrapped_cd = clock_drift.copy()
    valid_correlations = [
        c for c in bootstrapped_cd.correlations if not np.isnan(c.t_app[-1])]

    indices = np.random.choice(
        range(len(valid_correlations)),
        replace=True, size=len(valid_correlations))
    bootstrapped_cd.correlations = [valid_correlations[i] for i in indices]

    return bootstrapped_cd


def get_a_and_b_values(clock_drift):
    """
    Retrieve a and b values when the station was included in the inversion.
    """
    a_vals, b_vals, station_codes = [], [], []

    for station in clock_drift.stations:
        if station.needs_correction and (not all([station.a[-1] == 0, station.b[-1] == 0])):
            if any(not np.isnan(c.t_app[-1]) for c in clock_drift.get_correlations_of_station(
                    station.code)):
                a_vals.append(station.a[-1])
                b_vals.append(station.b[-1])
                station_codes.append(station.code)

    return a_vals, b_vals, station_codes


def perform_bootstrapping_and_solve(
        clock_drift, min_corr_params, days_apart,
        method='weighted_lstsq'):
    bootstrapped_cd = bootstrap_cd(clock_drift)

    with suppress_stdout():
        _ = bootstrapped_cd.filter_stations(*min_corr_params, days_apart)
        bootstrapped_cd.calculate_dt_ins()
        bootstrapped_cd.build_matrices()
        bootstrapped_cd.solve_eq(method=method)

    a_vals, b_vals, stations = get_a_and_b_values(bootstrapped_cd)

    if any(a > 1 for a in a_vals) or any(b > 1 for b in b_vals):
        return
    print(stations)

    for station_code, a, b in zip(stations, a_vals, b_vals):
        station = clock_drift.get_station(station_code)
        # This line updates the bootstrap_a attribute of the station object.
        # It uses the getattr function to safely access the attribute,
        # even if it doesn't exist initially.
        # If bootstrap_a doesn't exist, it defaults to an empty list [].
        station.bootstrap_a = getattr(station, 'bootstrap_a', []) + [a]
        station.bootstrap_b = getattr(station, 'bootstrap_b', []) + [b]
        station.bootstrap_inversion_method = method


def main(clock_drift, n_iterations=2):
    MIN_TOTAL_CORRELATIONS = 3
    MIN_CORRELATION_PERIODS = 2
    MIN_STATION_CONNECTIONS = 2
    DAYS_APART = 30

    min_corr_params = [MIN_TOTAL_CORRELATIONS,
                       MIN_CORRELATION_PERIODS, MIN_STATION_CONNECTIONS]

    for _ in range(n_iterations):
        perform_bootstrapping_and_solve(
            clock_drift, min_corr_params, DAYS_APART)


def restart_bootstrap(clock_drift):
    """Resets the bootstrap attributes for each station in the clock drift."""
    for station in clock_drift.stations:
        station.bootstrap_a = getattr(station, 'bootstrap_a', [])
        station.bootstrap_b = getattr(station, 'bootstrap_b', [])


with open("bootstrap/clockdrift_rebuttal_boots.obj", 'rb') as f:
    cd = pickle.load(f)

main(cd, n_iterations=2)

# %%
correlations = cd.get_correlations_of_station('O08')

# Using your provided commands
correlations_with_tapp = [c for c in correlations if not np.isnan(c.t_app[-1])]
avg_dates = [str(c.average_date)[:10] for c in correlations_with_tapp]

# Cluster dates
date_clusters = bootstrap_module.cluster_dates_by_range(
    avg_dates, days_apart=30)
print(date_clusters)

# %%


def meet_condition_total_correlations(correlations, min_number_of_total_correlations):
    """
    Check if a station meets the condition of having a minimum number of correlations with t_app.
    """
    no_correlations_with_t_app = sum(
        1 for c in correlations if not np.isnan(c.t_app[-1]))
    return no_correlations_with_t_app >= min_number_of_total_correlations


def meet_condition_correlation_periods(station, min_number_correlation_periods):
    """
    Check if a station meets the condition of having a minimum number of correlation periods.
    """
    return len(station.no_corr_per_avg_date) >= min_number_correlation_periods


def meet_condition_station_connections(correlations, min_number_of_stationconnections):
    """
    Check if a station meets the condition of having a minimum number of station connections.
    """
    connections_dictionary = {}
    for c in correlations:
        if not np.isnan(c.t_app[-1]):
            if station.code == c.station1_code:
                connections_dictionary[c.station2_code] = 1
            else:
                connections_dictionary[c.station1_code] = 1
    return len(connections_dictionary) >= min_number_of_stationconnections


def filter_stations(self,
                    min_number_of_total_correlations=3,
                    min_number_correlation_periods=2,
                    min_number_of_stationconnections=2,
                    days_apart=30):
    """
    Filter out stations that don't meet certain criteria.
    """
    stations_excluded_from_inversion = []

    for station in self.stations:
        if not station.needs_correction:
            continue

        if station in stations_excluded_from_inversion:
            continue

        # Recalculate the number of correlations per average date
        self.no_corr_per_avg_date(
            station=station, days_apart=days_apart, plot=False)

        correlations_of_station = self.get_correlations_of_station(
            station.code)

        # Check Conditions
        conditions = [
            meet_condition_total_correlations(
                correlations_of_station, min_number_of_total_correlations),
            meet_condition_correlation_periods(
                station, min_number_correlation_periods),
            meet_condition_station_connections(
                correlations_of_station, min_number_of_stationconnections)
        ]

        # If any condition is not met, set the t_app value to NaN for the station's correlations
        if not all(conditions):
            for c in correlations_of_station:
                c.t_app[-1] = np.nan
            stations_excluded_from_inversion.append(station)
            return False

    return True
