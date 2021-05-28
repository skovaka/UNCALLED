"""Simple example to read data from read until and send responses to server"""

import argparse
import functools
import logging
from multiprocessing.pool import ThreadPool
import time
import typing

import numpy

import read_until


def get_parser() -> argparse.ArgumentParser:
    """Build argument parser for example"""
    parser = argparse.ArgumentParser("Read until API demonstration..")
    parser.add_argument("--host", default="127.0.0.1", help="MinKNOW server host.")
    parser.add_argument(
        "--port", type=int, default=8000, help="MinKNOW gRPC server port."
    )
    parser.add_argument("--workers", default=1, type=int, help="worker threads.")
    parser.add_argument(
        "--analysis_delay",
        type=int,
        default=1,
        help="Period to wait before starting analysis.",
    )
    parser.add_argument(
        "--run_time", type=int, default=30, help="Period to run the analysis."
    )
    parser.add_argument(
        "--unblock_duration",
        type=float,
        default=0.1,
        help="Time (in seconds) to apply unblock voltage.",
    )
    parser.add_argument(
        "--one_chunk",
        default=False,
        action="store_true",
        help="Minimum read chunk size to receive.",
    )
    parser.add_argument(
        "--min_chunk_size",
        type=int,
        default=2000,
        help="Minimum read chunk size to receive. NOTE: this functionality "
        "is currently disabled; read chunks received will be unfiltered.",
    )
    parser.add_argument(
        "--debug",
        help="Print all debugging information",
        action="store_const",
        dest="log_level",
        const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        "--verbose",
        help="Print verbose messaging.",
        action="store_const",
        dest="log_level",
        const=logging.INFO,
    )
    return parser


def simple_analysis(
    client: read_until.ReadUntilClient,
    batch_size: int = 10,
    delay: float = 1,
    throttle: float = 0.1,
    unblock_duration: float = 0.1,
):
    """A simple demo analysis leveraging a `ReadUntilClient` to manage
    queuing and expiry of read data.

    :param client: an instance of a `ReadUntilClient` object.
    :param batch_size: number of reads to pull from `client` at a time.
    :param delay: number of seconds to wait before starting analysis.
    :param throttle: minimum interval between requests to `client`.
    :param unblock_duration: time in seconds to apply unblock voltage.

    """

    logger = logging.getLogger("Analysis")
    logger.warning(
        "Initialising simple analysis. "
        "This will likely not achieve anything useful. "
        "Enable --verbose or --debug logging to see more."
    )
    # we sleep a little simply to ensure the client has started initialised
    logger.info("Starting analysis of reads in %ss.", delay)
    time.sleep(delay)

    sample_count = 0

    while client.is_running:
        time_begin = time.time()
        # get the most recent read chunks from the client
        read_batch = client.get_read_chunks(batch_size=batch_size, last=True)
        for channel, read in read_batch:

            # convert the read data into a numpy array of correct type
            raw_data = numpy.frombuffer(read.raw_data, client.signal_dtype)
            read.raw_data = bytes("", "utf8")

            sample_count += len(raw_data)

            # make a decision that the read is good at we don't need more data?
            if (
                read.median_before > read.median
                and read.median_before - read.median > 60
            ):
                client.stop_receiving_read(channel, read.number)
            # we can also call the following for reads we don't like
            client.unblock_read(channel, read.number, duration=unblock_duration)

        # limit the rate at which we make requests
        time_end = time.time()
        if time_begin + throttle > time_end:
            time.sleep(throttle + time_begin - time_end)

    return sample_count


def run_workflow(
    client: read_until.ReadUntilClient,
    analysis_worker: typing.Callable[[], None],
    n_workers: int,
    run_time: float,
    runner_kwargs: typing.Optional[typing.Dict] = None,
):
    """Run an analysis function against a ReadUntilClient.

    :param client: `ReadUntilClient` instance.
    :param analysis worker: a function to process reads. It should exit in
        response to `client.is_running == False`.
    :param n_workers: number of incarnations of `analysis_worker` to run.
    :param run_time: time (in seconds) to run workflow.
    :param runner_kwargs: keyword arguments for `client.run()`.

    :returns: a list of results, on item per worker.

    """
    logger = logging.getLogger("Manager")

    if not runner_kwargs:
        runner_kwargs = {}

    results = []
    pool = ThreadPool(n_workers)
    logger.info("Creating %s workers", n_workers)
    try:
        # start the client
        client.run(**runner_kwargs)

        # start a pool of workers
        for _ in range(n_workers):
            results.append(pool.apply_async(analysis_worker))
        pool.close()

        # wait a bit before closing down
        time.sleep(run_time)
        logger.info("Sending reset")
        client.reset()
        pool.join()
    except KeyboardInterrupt:
        logger.info("Caught ctrl-c, terminating workflow.")
        client.reset()

    # collect results (if any)
    collected = []
    for result in results:
        try:
            res = result.get(3)
        except TimeoutError:
            logger.warning("Worker function did not exit successfully.")
            collected.append(None)
        except Exception:  # pylint: disable=broad-except
            logger.exception("Worker raise exception:")
        else:
            logger.info("Worker exited successfully.")
            collected.append(res)
    pool.terminate()
    return collected


def main(argv=None):
    """simple example main cli entrypoint"""
    args = get_parser().parse_args(argv)

    logging.basicConfig(
        format="[%(asctime)s - %(name)s] %(message)s",
        datefmt="%H:%M:%S",
        level=args.log_level,
    )

    read_until_client = read_until.ReadUntilClient(
        mk_host=args.host,
        mk_port=args.port,
        one_chunk=args.one_chunk,
        filter_strands=True,
    )

    analysis_worker = functools.partial(
        simple_analysis,
        read_until_client,
        delay=args.analysis_delay,
        unblock_duration=args.unblock_duration,
    )

    results = run_workflow(
        read_until_client,
        analysis_worker,
        args.workers,
        args.run_time,
        runner_kwargs={"min_chunk_size": args.min_chunk_size},
    )

    for idx, result in enumerate(results):
        logging.info("Worker %s received %s samples", idx + 1, result)
