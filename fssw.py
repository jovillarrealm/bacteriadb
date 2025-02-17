from collections import deque

def sliding_window(iterable, window_size):
    """
    Implements a fixed-size sliding window over an iterable.

    Args:
        iterable: The iterable to slide over (e.g., list, tuple, string, generator).
        window_size: The size of the sliding window.

    Yields:
        A tuple representing the current window.  If the iterable has fewer 
        elements than window_size before it is exhausted, then the final
        window will be padded with None values.
    """

    if not isinstance(window_size, int) or window_size <= 0:
        raise ValueError("Window size must be a positive integer.")

    window = deque()
    for item in iterable:
        window.append(item)
        if len(window) == window_size:
            yield tuple(window)  # Yield the current window as a tuple
            window.popleft() #Efficiently remove the oldest element
    
    # Handle the remaining elements if the iterable is shorter than window_size
    while window:
        if len(window) < window_size:
            window.append(None) #Pad with None values
        if len(window) == window_size:
            return tuple(window)
            #window.popleft()


if __name__ == "__main__":
    # Example Usage:

    data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    # Window size 3
    for window in sliding_window(data, 3):
        print(window)

    print("-" * 20)

    # Window size 2
    for window in sliding_window(data, 2):
        print(window)

    print("-" * 20)

    # Window size 4
    for window in sliding_window(data, 4):
        print(window)

    print("-" * 20)

    # Example with a string
    for window in sliding_window("abcdefgh", 3):
        print(window)

    print("-" * 20)

    #Example with iterable shorter than window size
    for window in sliding_window([1,2], 3):
        print(window)

    print("-" * 20)

    #Example with empty iterable
    for window in sliding_window([], 3):
        print(window)

    print("-" * 20)

    # Example with a generator (important for very large datasets)
    def my_generator(n):
        for i in range(n):
            yield i

    for window in sliding_window(my_generator(10), 3):
        print(window)