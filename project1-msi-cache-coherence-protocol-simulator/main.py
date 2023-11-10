import enum


# Cache state enum
class CacheStates(enum.Enum):
    invalid = "invalid"
    modified = "modified"
    shared = "shared"


# Read and Write hit and miss enum
class RWStatus(enum.Enum):
    RM = "read miss"
    RH = "read hit"
    WM = "write miss"
    WH = "write hit"


# Processor class: holds the cache of the processor
class Processor:
    def __init__(self, cache_size: int):
        self.cache = [None] * cache_size
        self.attached_bus = None
        self.state = CacheStates.invalid.value


# Bus class: holds a list of attached processors and the memory
class Bus:
    def __init__(self, memory_size: int):
        self.processors = []
        self.memory_size = memory_size
        self.memory = [None] * memory_size
        self.signal = None

    def add_processor(self, processor: Processor):
        processor.attached_bus = self
        self.processors.append(processor)

    def send_signal(self, signal: str):
        self.signal = signal

    def print_processor(self):
        return "pass"


# Returns the processor number for the parsed response
def get_processor_number(p: str):
    if response_segments[0] == "p1":
        return 1
    elif response_segments[0] == "p2":
        return 2
    elif response_segments[0] == "p3":
        return 3
    elif response_segments[0] == "p4":
        return 4
    else:
        return 0


def check_if_empty(memory: list):
    for location in memory:
        if location == None:
            return True
    return False


# Prints all the info
def print_info(bus: Bus):
    print("Memory: ", bus.memory)
    for i, processor in enumerate(bus.processors):
        print("Processor", i+1, "-- State:", processor.state, "-- Cache:", processor.cache)


# The provided processor snoops (checks) on the bus signal
def snoop(signal: str, processor: Processor):
    if signal == RWStatus.RM.value:  # read miss
        if processor.state == CacheStates.invalid.value:  # do nothing
            return
        elif processor.state == CacheStates.modified.value:  # Processor writes to memory
            bus.memory = processor.cache
            processor.state = CacheStates.shared.value
        elif processor.state == CacheStates.shared.value:  # do nothing
            return
    elif signal == RWStatus.RH.value:  # read hit
        if processor.state == CacheStates.invalid.value:  # do nothing
            return
        elif processor.state == CacheStates.modified.value:  # Processor writes to memory
            bus.memory = processor.cache
            processor.state = CacheStates.invalid.value
        elif processor.state == CacheStates.shared.value:  # Processor is now invalid
            processor.state = CacheStates.invalid.value
    elif signal == RWStatus.WM.value:  # write miss
        if processor.state == CacheStates.invalid.value:  # do nothing
            return
        elif processor.state == CacheStates.modified.value:  # Processor writes to memory
            bus.memory = processor.cache
            processor.state = CacheStates.shared.value
        elif processor.state == CacheStates.shared.value:  # change from shared to invalid
            processor.state = CacheStates.invalid.value
    elif signal == RWStatus.WH.value:  # write hit
        if processor.state == CacheStates.invalid.value:  # do nothing
            return
        elif processor.state == CacheStates.modified.value:  # Processor writes to memory
            bus.memory = processor.cache
            processor.state = CacheStates.invalid.value
        elif processor.state == CacheStates.shared.value:  # Processor no longer valid
            processor.state = CacheStates.invalid.value


# Reading operatin
def read(bus: Bus, response_segments: list):
    # Setting up some variables to use
    processor_number = get_processor_number(response_segments[0])  # store p1-p4 as int 1-4
    processor = bus.processors[processor_number - 1]  # store the processor based on processor_number
    memory_location = int(response_segments[2]) - 1  # store the provided memory location (not needed since whole block is returned, and we are dealing with only 1 block)

    # Do stuff...
    if not processor_number == 0:  # sanity check
        if processor.state == CacheStates.invalid.value:  # invalid state
            bus.send_signal(RWStatus.RM.value)  # puts read miss on bus
            for i, p in enumerate(bus.processors):  # other processors snoop
                if not i == processor_number - 1:  # skip current processor
                    snoop(RWStatus.RM.value, p)
            processor.cache = bus.memory  # transfer data from memory to cache
            processor.state = CacheStates.shared.value  # change processor state
        elif processor.state == CacheStates.modified.value:  # modified state: do nothing, already has correct data
            pass
        elif processor.state == CacheStates.shared.value:  # shared state: do nothing, already has correct data
            pass

    # Print the results
    print_info(bus)


# Writing operation
def write(bus: Bus, response_segments: list):
    # Setting up some variables to use
    processor_number = get_processor_number(response_segments[0])
    processor = bus.processors[processor_number - 1]
    memory_location = int(response_segments[2]) - 1
    value = int(response_segments[3])

    # Do stuff...
    if not processor_number == 0:  # sanity check
        if processor.state == CacheStates.invalid.value:  # invalid state
            bus.send_signal(RWStatus.WM.value)  # puts write miss on bus
            for i, p in enumerate(bus.processors):  # other processors snoop
                if not i == processor_number - 1:  # skip current processor
                    snoop(RWStatus.WM.value, p)
            processor.cache = bus.memory  # transfer data from memory to cache
            processor.cache[memory_location] = value  # write value to cache
            processor.state = CacheStates.shared.value  # change processor state
        elif processor.state == CacheStates.modified.value:  # modified state: just write value to cache, already has correct data
            bus.send_signal(RWStatus.WH.value)  # puts write hit on bus
            for i, p in enumerate(bus.processors):  # other processors snoop
                if not i == processor_number - 1:  # skip current processor
                    snoop(RWStatus.WH.value, p)
            processor.cache[memory_location] = value  # write value to cache
            processor.state = CacheStates.modified.value  # change processor state
        elif processor.state == CacheStates.shared.value:  # shared state
            bus.send_signal(RWStatus.WH.value)  # puts write hit on bus
            for i, p in enumerate(bus.processors):  # other processors snoop
                if not i == processor_number - 1:  # skip current processor
                    snoop(RWStatus.WH.value, p)
            processor.cache[memory_location] = value  # write value to cache
            processor.state = CacheStates.modified.value  # change processor state

    # Print the results
    print_info(bus)


# Main "function"
if __name__ == "__main__":
    print("Initializing processors and memory...")

    # Setting up processors and bus; attaching processors to bus
    p1 = Processor(4)
    p2 = Processor(4)
    p3 = Processor(4)
    p4 = Processor(4)
    bus = Bus(4)
    bus.add_processor(p1)
    bus.add_processor(p2)
    bus.add_processor(p3)
    bus.add_processor(p4)

    # Putting some initial values in memory
    bus.memory[0] = 10
    bus.memory[1] = 20
    bus.memory[2] = 30
    bus.memory[3] = 40

    # Print the initial states
    print_info(bus)

    # Declaring response variable
    quit = False

    # Main loop
    while not quit:
        # Show some info and wait for response
        response = input(
            "\n------------------------------\nReading:\t[p1, p2, p3, p4] [read, load] [memory location]\nWriting:\t[p1, p2, p3, p4] [write, store] [memory_location] [value]\nmemory_location = 1, 2, 3, or 4\nQuitting:\tQ or q\n------------------------------\nWhat do you want to do? ")
        print("\n")

        if response == "q" or response == "Q":  # Quit
            exit()
        else:
            # Split response into usable segments
            response_segments = response.split(' ')

            # Do some brief sanity checks
            if len(response_segments) >= 3 and len(response_segments) <= 4:
                if not (response_segments[0] == "p1" or response_segments[0] == "p2" or response_segments[0] == "p3" or response_segments[0] == "p4"):
                    print("First argument is not valid")
                    continue
                if not (response_segments[1] == "read" or response_segments[1] == "load" or response_segments[1] == "write" or response_segments[1] == "store"):
                    print("Second argument is not valid")
                    continue
                if not response_segments[2].isdigit():
                    print("Third argument is not an int")
                    continue
                if len(response_segments) == 4 and (response_segments[1] == "read" or response_segments[1] == "load"):
                    if not response_segments[3].isdigit():
                        print("Fourth argument is not an int")
                        continue
            else:
                print("Invalid number of arguments")
                continue

            # Handling reading and writing operations
            if response_segments[1] == "read" or response_segments[1] == "load":  # Do a read operation
                read(bus, response_segments)
            elif response_segments[1] == "write" or response_segments[1] == "store":  # Do a write operation
                write(bus, response_segments)
