ifndef inputs
$(error "Project has no defined $${inputs}.")
endif
ifndef outputs
$(error "Project has no defined $${outputs}.")
endif
ifndef root
$(error "The $${root} path is not defined.")
endif
ifndef tooldir
$(error "The $${tooldir} path is not defined.")
endif
ifndef remotedir
$(error "The $${remotedir} path is not defined.")
endif
tools += dispatcher.py worker.py
userscripts = prepare.py extract.py process.py
external += ${tools} proj.mk
remotedir := ${remotedir}/$(PWD:$(wildcard ${root})/%=%)

.SECONDEXPANSION:
.PRECIOUS: $(addprefix results_%/,${outputs}) results_%/results.p RUN/%_job.log

local:
	@${MAKE} --no-print-directory $(addprefix results_$@/,${outputs})

$(addprefix results_%/,${outputs}): results_%/results.p process.py
	cd results_$* && python ../process.py

results_%/results.p: RUN/%_job.log extract.py | ${external}
ifneq ("$(wildcard RUN/*.start RUN/*.running.*)", "")
	$(warning "Some jobs are still running.")
endif
	python extract.py
	mkdir -p results_$* && mv RUN/results.p $@

RUN/%_job.log: prepare.py ${inputs}
	@${MAKE} --no-print-directory prepare
	@${MAKE} --no-print-directory run_$*
	@$(if $(subst local,,$*), @${MAKE} --no-print-directory print_error)

print_error:
	$(error "Wait till the job finishes, then run make again.")

${tools} proj.mk:
	@mkdir -p $(dir $@)
	@rsync -ai ${tooldir}/$@ $(dir $@)

run_local:
	python worker.py RUN 1 >RUN/local_job.log

run_%:
	bash ~/bin/submit.sh $*.job.sh
	@sleep 1  # some submitters print asynchronously
	
prepare: | ${external}
ifneq ("$(wildcard RUN)", "")
	$(error "There is a previous RUN, run make cleanrun to overwrite.")
endif
	${prepare_env} python prepare.py

update:
	@echo "Updating tools..."
ifneq ("$(wildcard ${external})", "")
	@rsync -aR $(wildcard ${external}) .oldtools/
endif
	@${MAKE} --no-print-directory -B external

restore:
	@echo "Restoring tools..."
	@rsync -ia .oldtools/* ./

external: ${external}

remote_%: upload_$$(firstword $$(subst _, , %)) 
ifdef OFFLINE
	@echo "Skipping download."
else
	$(eval remote := $(firstword $(subst _, ,$*)))
	@echo "Connecting to ${remote}..."
	@ssh ${remote} \
		"cd ${remotedir} && \
		${prepare_env_remote} make results_$*/results.p"
	@echo "Downloading results from ${remote}..."
	@rsync -ia ${remote}:${remotedir}/results_$*/results.p results_$*/
endif
	@${MAKE} --no-print-directory $(addprefix results_$*/,${outputs})

upload_%: ${external}
ifdef OFFLINE
	@echo "Skipping upload."
else
	@echo "Uploading to $*..."
	@ssh $* "mkdir -p ${remotedir}"
	@rsync -ia \
		--exclude=*.pyc --exclude=RUN $(addprefix --exclude=,${excluded}) \
		--include=$*_*.job.sh --exclude=*_*.job.sh \
		--exclude=results_* \
		${PWD}/* $*:${remotedir}/
endif

submit_%:
	$(eval remote := $(firstword $(subst _, ,$*)))
	@echo "Connecting to ${remote}..."
	@ssh ${remote} "cd ${remotedir} && make run_$*"

archive_%:
	@${MAKE} --no-print-directory results_$*/$(notdir ${PWD}).tar.gz

results_%/$(notdir ${PWD}).tar.gz: $$(addprefix results_%/,${outputs})
	@echo "Creating archive $@..."
	@tar -zc ${tools} ${userscripts} $(addprefix results_$*/,${outputs}) ${inputs} \
		${external} Makefile >$@

clean:
ifneq ("$(wildcard *.pyc)", "")
	rm *.pyc
endif

cleanrun:
ifneq ("$(wildcard RUN)", "")
	rm -r RUN
endif

distclean: clean cleanrun
ifneq ("$(wildcard ${tools} ${excluded} results_*/*)", "")
	rm -r $(wildcard ${tools} ${excluded} results_*)
endif
	-rm -rf .oldtools

monitor_%:
	@ssh $* qmy

check: numoftasks

numoftasks:
	@echo "Number of initialized tasks: $(shell ls -d RUN/*.start 2>/dev/null | wc -l)"
	@echo "Number of running tasks: $(shell ls -d RUN/*.running.* 2>/dev/null | wc -l)"
	@echo "Number of finished tasks: $(shell ls -d RUN/*.done 2>/dev/null | wc -l)"

check_%:
	@echo "Connecting to $*..."
	@ssh $* "cd ${remotedir} && make check"

cleanrun_%:
	@echo "Connecting to $*..."
	@ssh $* "cd ${remotedir} && make cleanrun"

distclean_%:
	@echo "Connecting to $*..."
	@ssh $* "cd ${remotedir} && make distclean"
