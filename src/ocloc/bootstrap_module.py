from datetime import datetime, timedelta


def convert_strings_to_dates(date_strings):
    """Convert a list of date strings to a sorted list of datetime.date objects."""
    return sorted([datetime.strptime(date_str, "%Y-%m-%d").date() for date_str in date_strings])


def get_average_date(cluster_dates):
    """Get the average date for a list of dates."""
    # Calculate average date by averaging the timestamps
    avg_timestamp = sum([date.toordinal()
                        for date in cluster_dates]) / len(cluster_dates)
    avg_date = datetime.fromordinal(int(avg_timestamp))
    return avg_date.date()


def cluster_dates_by_range(dates, days_apart=10):
    """
    Group dates that are within a certain range of each other.

    Parameters:
    - dates (list of str): The dates to be clustered.
    - days_apart (int): The maximum number of days between dates in a cluster.

    Returns:
    - dictionary: Keys are the average dates and values are the number of dates in each cluster.
    """
    sorted_dates = convert_strings_to_dates(dates)
    clusters = []
    current_cluster = [sorted_dates[0]]

    for date in sorted_dates[1:]:
        if date - current_cluster[-1] <= timedelta(days=days_apart):
            current_cluster.append(date)
        else:
            clusters.append(current_cluster)
            current_cluster = [date]

    clusters.append(current_cluster)

    average_dates_dict = {}
    for cluster in clusters:
        avg_date = get_average_date(cluster)
        average_dates_dict[avg_date.strftime("%Y-%m-%d")] = len(cluster)

    return average_dates_dict
