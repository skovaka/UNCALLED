import traceback
from threading import Event, Thread
import sys
import time
import numpy as np
import logging

try:
    import read_until
    ru_loaded = True
except ImportError:
    ru_loaded = False

if ru_loaded:

    MK_ERROR_STATUS = 0
    MK_READY = 1
    MK_STARTING = 2
    MK_PROCESSING = 3
    MK_FINISHING = 4

    class Client(read_until.ReadUntilClient):
        def __init__(self, 
                     mk_host=8000, 
                     mk_port="127.0.0.1", 
                     chunk_len=4000, 
                     scan_thresh=0.99):

            logging.basicConfig(
                format='[%(asctime)s - %(name)s] %(message)s',
                datefmt='%H:%M:%S', level=logging.INFO)

            self.chmon_running = Event()
            self.chmon_thread = None

            read_until.ReadUntilClient.__init__(
                    self,
                    mk_host=mk_host,  
                    mk_port=mk_port,  
                    one_chunk=False,    
                    filter_strands=True,
                    prefilter_classes={'strand', 'strand1', 'strand2'}
            )


            self.in_scan = True 
            self.scan_thresh = scan_thresh
            self.chunk_len = chunk_len

            self.ch_mux = np.zeros(512, dtype=int)
            self.mux_counts = np.zeros(5, dtype=float)
            self.mux_counts[0] = len(self.ch_mux)

            self.logger.setLevel(logging.WARNING)
            self.unc_log = logging.getLogger('UNCALLED')


        def run(self, steady_wait=10, scan_wait=5, refresh=0.5, **kwargs):

            self._wait_for_start(steady_wait, refresh)

            self._start_chmon()

            self._set_chunk_len(scan_wait)

            read_until.ReadUntilClient.run(self, *kwargs)

            self.start_time = time.time()

        def reset(self, timeout=5):
            read_until.ReadUntilClient.reset(self, timeout)

            if self.chmon_thread is not None:
                self.chmon_running.clear()
                self.chmon_thread.join()
                if self.chmon_thread.is_alive():
                    self.unc_log.error("Channel monitor didn't exit right")
            self.chmon_thread = None

        def get_runtime(self):
            return time.time() - self.start_time

        def should_eject(self):
            return not self.in_scan
        
        def log(self, msg):
            self.unc_log.info(msg)

        def _start_chmon(self):
            self.start_time = time.time()
            self.last_scan = time.time()
            self.in_scan = True

            self.chmon_running.set()
            self.chmon_thread = Thread(
                target=self._run_chmon,
                name="channel_monitor"
            )
            self.chmon_thread.start()

        def _scan_update(self):
            #ps = self.in_scan

            m = np.argmax(self.mux_counts)
            self.in_scan = (m != 0 and 
                            (self.mux_counts[m] / (len(self.ch_mux)-self.mux_counts[0])
                                > self.scan_thresh))

            if self.in_scan:
                self.last_scan = time.time()

            #if self.in_scan != ps:
            #    sys.stderr.write("%.2f in_scan changed to %d\n" % 
            #                     (self.get_runtime(), self.in_scan))

        def _get_status(self):
            return self.connection.acquisition.current_status().status

        def _wait_for_start(self, steady=10, refresh=0.5):
            if self._get_status() == MK_PROCESSING:
                self.log("Run already in progress")
                return True

            processing = False
            proc_time = None
            
            self.log("Waiting for run to start")

            while True:
                if self._get_status() == MK_PROCESSING:
                    if proc_time == None:
                        proc_time = time.time()

                        if not processing:
                            self.log("Waiting for steady state")
                            processing = True

                    elif time.time() - proc_time >= steady:
                        return True

                elif proc_time != None:
                    proc_time = None

                time.sleep(refresh)

            self.last_scan = time.time()

        def _set_chunk_len(self, scan_wait=5):

            #Don't want to restart acquisition during mux scan
            scan_gap = time.time() - self.last_scan
            if self.in_scan or scan_gap < scan_wait:
                self.log("Checking for initial mux scan")
                while self.in_scan or scan_gap < scan_wait:
                    time.sleep(max(0.1, scan_wait - scan_gap))
                    scan_gap = time.time() - self.last_scan

            chunk_sec = self.chunk_len/4000.0

            client = self.connection.analysis_configuration
            config = client.get_analysis_configuration()

            prev_sec = config.read_detection.break_reads_after_seconds.value

            self.log("Setting chunk length to %.2f seconds" % chunk_sec)
            if prev_sec != chunk_sec:

                self.connection.acquisition.stop(wait_until_ready=True, keep_power_on=True)

                config.read_detection.break_reads_after_seconds.value=chunk_sec
                client.set_analysis_configuration(config)

                self.connection.acquisition.start(wait_until_processing=True)

                self._start_chmon()

        def _update_muxs(self, channels):
            for channels in channels:
                if not self.chmon_running.is_set():
                    break

                for ch in channels.channel_states:
                    nmx = int(ch.config.well)
                    omx = self.ch_mux[ch.channel-1]
                    if omx != nmx:
                        self.mux_counts[omx] -= 1
                        self.mux_counts[nmx] += 1
                        self.ch_mux[ch.channel-1] = nmx

                self._scan_update()

                if self._get_status() != MK_PROCESSING:
                    self.chmon_running.clear()
                    self.running.clear()
                    break


        def _run_chmon(self, **kwargs):
            channels = self.connection.data.get_channel_states(
                first_channel = 1,
                last_channel = 512,
                use_channel_states_ids = False
            )

            try:
                self._update_muxs(channels)
            except Exception as e:
                sys.unc_log.error(traceback.format_exc())

            # signal to the server that we are done with the stream.
            channels.cancel()

