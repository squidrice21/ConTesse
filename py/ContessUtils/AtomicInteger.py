# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

# ADOBE CONFIDENTIAL
import threading

class AtomicInteger:
    def __init__(self):
        self._lock = threading.Lock()
        self._value = 0

    def increment_and_get(self):
        with self._lock:
            self._value += 1
            return self._value

    def decrement_and_get(self):
        with self._lock:
            self._value -= 1
            return self._value

    def get(self):
        with self._lock:
            return self._value