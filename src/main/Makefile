BUILD = $(shell find `pwd` -mindepth 1 -maxdepth 1 -type d)

lib:
	$(foreach folder, $(BUILD), mkdir -p $(folder)/build;)
	$(foreach folder, $(BUILD), cd $(folder); $(MAKE);)

.PHONY: clean
clean:
	$(foreach folder, $(BUILD), cd $(folder); $(MAKE) clean;)