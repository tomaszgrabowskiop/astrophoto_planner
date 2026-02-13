import time
import json

class SmartCache:
    def __init__(self, expiration_time=3600):
        self.cache = {}
        self.expiration_time = expiration_time

    def get(self, key):
        if key in self.cache:
            value, timestamp = self.cache[key]
            if time.time() - timestamp < self.expiration_time:
                return value
            else:
                del self.cache[key]
                return None
        return None

    def set(self, key, value):
        self.cache[key] = (value, time.time())


def smart_cache(func):
    cache = SmartCache()

    def wrapper(*args):
        key = json.dumps(args)  # Create a cache key based on function arguments
        cached_result = cache.get(key)
        if cached_result is not None:
            return cached_result
        
        result = func(*args)
        cache.set(key, result)
        return result

    return wrapper

# Example usage:
@smart_cache

def expensive_computation(a, b):
    time.sleep(2)  # Simulate a time-consuming computation
    return a + b

# Example of using the function with caching
print(expensive_computation(1, 2))  # This will take time to compute
print(expensive_computation(1, 2))  # This will return immediately from cache
